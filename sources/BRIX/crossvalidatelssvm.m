function cost = crossvalidatelssvm(model,Y, L, omega, estfct,combinefct)

nb_data = size(Y,1);
d = size(model.xtrain,2);
% LS-SVMlab
eval('model = initlssvm(model{:});',' ');
model.status = 'changed';

eval('L;','L=min(round(sqrt(size(model.xfull,1))),10);');
eval('estfct;','estfct=''mse'';');
eval('combinefct;','combinefct=''mean'';');

% Y is raw data, non preprocessed
py = Y;
[~,Y] = postlssvm(model,[],Y);

gams = model.gamcsa; try sig2s = model.kernel_parscsa; catch, sig2s = [];end

%initialize: no incremental  memory allocation
costs = zeros(L,length(gams));
block_size = floor(nb_data/L);

% check whether there are more than one gamma or sigma
for j =1:numel(gams)
    if strcmp(model.kernel_type,'RBF_kernel') || strcmp(model.kernel_type,'RBF4_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',sig2s(j));
    elseif strcmp(model.kernel_type,'lin_kernel')
        model = changelssvm(model,'gam',gams(j));
    elseif strcmp(model.kernel_type,'poly_kernel')
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',[sig2s(1,j);sig2s(2,j)]);
    else
        model = changelssvm(changelssvm(model,'gam',gams(j)),'kernel_pars',[sig2s(1,j);sig2s(2,j);sig2s(3,j)]);
    end
    
    
    % calculate matrix for LS-SVM once for the entire data
    S = ones(nb_data,1);
    Inb = eye(nb_data);
    K = kernel_matrix2(omega,model.kernel_type,model.kernel_pars,d);
    Atot = K+Inb./model.gam;

    %profile -memory on
    profile on
    
        A = [K+Inb./model.gam S; S' 0];
        C = pinv(A);
        alpha = C*[py;0];
        alpha = alpha(1:nb_data);
        clear A
    
    % clear C
    clear Atot
    block_size = (nb_data+1)/L;
    disp(block_size); pause;
   
    for l = 1:L,
        if l==L,
            validation = 1:block_size-1;
        else
            validation = 1:block_size;
        end
        
        BRI_in = [L block_size d model.gam model.kernel_pars l-1];
        Ckk = BRI(model.xtrain,BRI_in);
        % disp(Ckk); pause;
        Ckk = Ckk(validation,validation);
        % disp(Ckk); pause;
        X = model.xtrain;
        whos C
        whos X
        whos Ckk
        
        profile viewer
        pause;
        
        try % faster
            Rkk = chol(Ckk+eps);
            betak = Rkk\(Rkk'\alpha(validation));
        catch
            betak = Ckk\alpha(validation);
        end
        % latent outputs for validation
        yh = py(validation) - betak;
        [~,yh] = postlssvm(model,[],yh);
        if ~(model.type(1)=='c')
            costs(l,j) = feval(estfct,yh - Y(validation,:));
        else
            costs(l,j) = feval(estfct,Y(validation,:),sign(yh));
        end
    end
end
cost = feval(combinefct, costs);
