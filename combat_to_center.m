
%data = randn(1000,10);
%batch = [1 1 1 1 1 2 2 2 2 2];
%mod = [1 2 1 2 1 2 1 2 1 2]';



function bayesdata = combat_to_center(dat, batch, mod, parametric,center)
    [sds] = std(dat')'; %得到每一个脑区的标准差，假设n个脑区，这size(sds) = [n,1]
    wh = find(sds==0); %检验是否存在标准差为0的一列，若为零，则说明存在一列数是常数，程序会有错误
    [ns,ms] = size(wh);
    if ns>0
        error('Error. There are rows with constant values across samples. Remove these rows and rerun ComBat.')
    end
    batchmod = categorical(batch);
    batchmod = dummyvar({batchmod});%这两部是将其变为分类数据，并转化为哑变量，即类别1变为1,0，类别2变为0，1等等
	n_batch = size(batchmod,2);%n_batch是batch的数目
	levels = unique(batch);%levels里记录了batch的原始编号，用于得到每个batch的位置
	fprintf('[combat] Found %d batches\n', n_batch);

	batches = cell(0);%batches是一个cell数组，用于记录每个batch的位置，一共有几个batch就有几个cell
	for i=1:n_batch
		batches{i}=find(batch == levels(i));
	end
	n_batches = cellfun(@length,batches);%n_batches记录了每个batch内被试的个数
	n_array = sum(n_batches);%所有被试的个数

	% Creating design matrix and removing intercept:
	design = [batchmod mod];
	intercept = ones(1,n_array)';
	wh = cellfun(@(x) isequal(x,intercept),num2cell(design,1));%cellfun(func,C)是对元胞数组C中每个元胞应用函数func
	bad = find(wh==1);
	design(:,bad)=[];%假设design中存在截距项，即全是1的值，将其去掉


	fprintf('[combat] Adjusting for %d covariate(s) of covariate level(s)\n',size(design,2)-size(batchmod,2))
	% Check if the design is confounded
	if rank(design)<size(design,2) %rank是矩阵的秩，矩阵的秩是指矩阵中线性无关的行或列的最大数目
		nn = size(design,2);
	    if nn==(n_batch+1) 
	      error('Error. The covariate is confounded with batch. Remove the covariate and rerun ComBat.')
	    end
	    if nn>(n_batch+1)
	      temp = design(:,(n_batch+1):nn);
	      if rank(temp) < size(temp,2)
	        error('Error. The covariates are confounded. Please remove one or more of the covariates so the design is not confounded.')
	      else 
	        error('Error. At least one covariate is confounded with batch. Please remove confounded covariates and rerun ComBat.')
	      end
	    end
	 end


	fprintf('[combat] Standardizing Data across features\n')
	B_hat = inv(design'*design)*design'*dat';%对dat = design*B_hat + err用最小二乘法估计B_hat
	%Standarization Model
	grand_mean = (n_batches/n_array)*B_hat(1:n_batch,:);%等价于每个batch对应beta的加权平均，权重为每个batch对应的个体数目
	var_pooled = ((dat-(design*B_hat)').^2)*repmat(1/n_array,n_array,1); %等价于var_pooled = mean(((dat-(design*B_hat)').^2),2);
	stand_mean = grand_mean'*repmat(1,1,n_array);%此时stand_mean应该对应的是公式中Ag的估计值，Ag is the overal gene expression
	% Making sure pooled variances are not zero:
	wh = find(var_pooled==0);
	var_pooled_notzero = var_pooled;
	var_pooled_notzero(wh) = [];
	var_pooled(wh) = median(var_pooled_notzero);

	if not(isempty(design))
		tmp = design;
		tmp(:,1:n_batch) = 0;
		stand_mean = stand_mean+(tmp*B_hat)';%此时stand_mean应该是公式中Ag+Xbeta的值
	end	
	s_data = (dat-stand_mean)./(sqrt(var_pooled)*repmat(1,1,n_array));

	%Get regression batch effect parameters
	fprintf('[combat] Fitting L/S model and finding priors\n') % L/S代表了Location and scale
	batch_design = design(:,1:n_batch);
	gamma_hat = inv(batch_design'*batch_design)*batch_design'*s_data';
	delta_hat = [];
	for i=1:n_batch
		indices = batches{i};
		delta_hat = [delta_hat; var(s_data(:,indices)')];
	end

	%Find parametric priors:
	gamma_bar = mean(gamma_hat');
	t2 = var(gamma_hat');
	delta_hat_cell = num2cell(delta_hat,2);
	a_prior=[]; b_prior=[];
	for i=1:n_batch
		a_prior=[a_prior aprior(delta_hat_cell{i})];
		b_prior=[b_prior bprior(delta_hat_cell{i})];
	end

	
	if parametric
        fprintf('[combat] Finding parametric adjustments\n')
        gamma_star =[]; delta_star=[];
        for i=1:n_batch
            indices = batches{i};
            temp = itSol(s_data(:,indices),gamma_hat(i,:),delta_hat(i,:),gamma_bar(i),t2(i),a_prior(i),b_prior(i), 0.001);
            gamma_star = [gamma_star; temp(1,:)];
            delta_star = [delta_star; temp(2,:)];
        end
    end
	    
    if (1-parametric)
        gamma_star =[]; delta_star=[];
        fprintf('[combat] Finding non-parametric adjustments\n')
        for i=1:n_batch
            indices = batches{i};
            temp = inteprior(s_data(:,indices),gamma_hat(i,:),delta_hat(i,:));
            gamma_star = [gamma_star; temp(1,:)];
            delta_star = [delta_star; temp(2,:)];
        end
    end
	    
	fprintf('[combat] Adjusting the Data\n')
	bayesdata = s_data;
	j = 1;
	for i=1:n_batch
		indices = batches{i};
		bayesdata(:,indices) = (bayesdata(:,indices)-(batch_design(indices,:)*gamma_star)')./(sqrt(delta_star(j,:))'*ones(1,n_batches(i)));
		j = j+1;
    end
    
    for i=1:n_batch
		indices = batches{i};
		bayesdata(:,indices) = (bayesdata(:,indices)).*(sqrt(delta_star(center,:))'*ones(1,n_batches(i))) + gamma_star(center,:)'*ones(1,n_batches(i));
    end
    
	bayesdata = (bayesdata.*(sqrt(var_pooled)*ones(1,n_array)))+stand_mean;

end

