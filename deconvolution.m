%% jMF2D deconvolution
function[F,Z]=deconvolution(X,Y,L2,alpha,beta,gamma,loss,max_iter, ...
    seed,initial_B,initial_F, sym, normbyspot)
    %% initial
    s = size(X,2);
    g = size(X,1); 
    k = size(Y,2); 
    if seed~=0
        rng(seed);
    end
    if initial_B==0
        B = rand(g,k);
    elseif initial_B==1
        B = Y;
    else
        B = Y +rand(g,k);
    end
    if initial_F==0
        F = rand(k,s);
    else
        F = B \ X;
    end
    E = rand(k,s);
    Z = ones(k,k);
    C = ones(k,k);
    T1 = zeros(k,k);
    T2 = zeros(k,s);
    delta1 = 1; 
    delta2 = 1;
    
    %% Update parameters iteratively
    converge = 0; iter = 0;
    tic;
    while ~converge && iter<max_iter
        % Update Z
        Zk = Z;
        Z = Z - diag(diag(Z));
        Z = Z.*((Y'*Y+delta1*C-T1)./(Y'*Y*Z+delta1*Z));
        Z = Z./repmat(sum(Z,2),1,size(Z,2));
        Z = Z - diag(diag(Z)-1);
        if sym==0
            Z = (Z + Z') / 2;
        else
            Z = triu(Z);
            Z = Z + Z' - eye(k);
        end
        tz = norm(Z-Zk,'fro');
        % cal L1
        L1 = diag(sum(Z))-Z;
        % Update C
        Ck = C;
        Crow = sqrt(sum(abs(conj(C').*C')));
        D1 = diag(1/2*Crow);
        C = C.*((delta1*Z+T1)./(alpha*D1*C+delta1*C));
        tc = norm(C-Ck,'fro');
        % Update T1
        T1 = T1 + delta1*(Z-C);
        % Update B
        Bk = B;
        B = B.*((X*F')./(B*F*F'+(beta/2)*(B*L1+B*L1')));
        tb = norm(B-Bk,'fro');
        B = max(B,0);
        B(find(isnan(B)==1)) = 0;
        % Update F
        Fk = F;
        F = F.*((B'*X+delta2*E)./(B'*B*F+delta2*F+T2));
        if normbyspot==1
            F = F./repmat(sum(F,1),size(F,1),1);
        end
        F = real(F);
        tf = norm(F-Fk,'fro');
        F = max(F,0);
        % Update E
        Ek = E;
        Erow = sqrt(sum(abs(conj(E').*E')));
        D2 = diag(1/2*Erow);
        E = E.*((delta2*F+T2)./((gamma/2)*(E*L2+E*L2')+alpha*D2*E+delta2*E));
        E = max(E,0);
        te = norm(E-Ek,'fro');
        % Update T2
        T2 = T2 + delta2*(F-E);
        % cal loss
        tol_loss = max([norm(B-Bk,'fro'),tf]) / max([norm(B,'fro'),norm(F,'fro')]);
        if tol_loss < loss
            converge = 1;
        else
            iter = iter + 1;
        end
        disp([' iter ' num2str(iter) ' loss ' num2str(tol_loss)]);
    end
    disp(['Runtime (s):' num2str(toc)]);
end