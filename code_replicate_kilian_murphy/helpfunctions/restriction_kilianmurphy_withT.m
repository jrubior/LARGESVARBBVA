function [St_old] = restriction_kilianmurphy_withT(Bdraw, Sigmadraw, Qdraw,fo_inv,fo_str2irfs,info)
% checking sign restrictions


    %BSigmaQnew = [vec(Bdraw);vec(Sigmadraw);vec(Qdraw)];

    %structpara   = fo_inv(BSigmaQnew);

    %irfpara      = fo_str2irfs(structpara);




     Rr=chol(Sigmadraw)'*Qdraw;


    Tt=zeros(info.nshocks,info.nvar);
    for j=1:info.nshocks

        switch j



            case 1

                for i=1:info.nvar
                    if all(info.Ss{j,1}*Rr(:,i)>0)

                        Tt(j,i)=1;
                    elseif all(info.Ss{j,1}*Rr(:,i)<0)
                        Tt(j,i)=-1;
                    end
                end


            case 2

                for i=1:info.nvar
                    if all(info.Ss{j,1}*Rr(:,i)>0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=1;
                    elseif all(info.Ss{j,1}*Rr(:,i)<0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=-1;
                    end
                end

            case 3

                for i=1:info.nvar
                    if all(info.Ss{j,1}*Rr(:,i)>0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=1;
                    elseif all(info.Ss{j,1}*Rr(:,i)<0) && ((Rr(1,i)/Rr(3,i))<0.0258)
                        Tt(j,i)=-1;
                    end

                end

        end
    end

    checkA1 = any(sum(Tt ~= 0, 1) > 1);

    if checkA1==1
        disp('PROBLEM with Tt')
        return;
    end


    %if rank(Tt)<nshocks
    
    if any(sum(Tt == 0, 2)==size(Tt,2))
        label=0;
    else


        Rstar = zeros(info.nvar,info.nvar);
        S = cell(info.nshocks,1);

        isamples = zeros(info.nshocks,1);

        % For each row j
        for j = 1:info.nshocks
      
            S{j} = find(Tt(j,:) == 1 | Tt(j,:) == -1);



            % Get a random index from S_j
            if numel(S{j})>1
                randomIndex = randsample(S{j},1);
            else
                randomIndex = S{j};
            end

            isamples(j,1) = randomIndex;



            if Tt(j,isamples(j))==1
                Rstar(:,j) = Rr(:,isamples(j,1));
            elseif Tt(j,isamples(j))==-1
                Rstar(:,j) = -Rr(:,isamples(j,1));
            end

        end

        % Loop through j = m+1, ..., n
        i_values=zeros(info.nvar-info.nshocks,1);
        for j = info.nshocks+1:info.nvar
            % Create S_j = {1,...,n} \ {i_1,...,i_j}
            % We need to exclude i_1 through i_(j-1) from the set {1,...,n}
           
            S_j = setdiff(1:info.nvar, isamples(1:j-1,1));

            % Sample an element uniformly from S_j
            if numel(S_j)>1


                random_index = randsample(S_j,1);

            else
                random_index = S_j;
            end
            sampled_element = random_index;

            % Store the sampled element as i_j
            i_values(j) = sampled_element;

            if rand(1)>0.5
                Rstar(:,j) = Rr(:,i_values(j));
            else
                Rstar(:,j) = -Rr(:,i_values(j));
            end

        end

        Qdraw=(chol(Sigmadraw)')\Rstar;

        L0old=Rstar;


            %check dynamic sign restrictions
                    % dynamic impulse responses
            hSigmadraw = chol(Sigmadraw);
            A0         = hSigmadraw\Qdraw;
            Aplus      = Bdraw*A0;
            for l=1:info.nlag-1
                info.A{l} = Aplus((l-1)*info.nvar+info.nex+1:l*info.nvar+info.nex,1:end);
                info.bigF((l-1)*info.nvar+1:l*info.nvar,1:info.nvar)=info.A{l}/A0;
            end
            A{info.nlag} = Aplus((info.nlag-1)*info.nvar+1:info.nlag*info.nvar,1:end);
            info.bigF((info.nlag-1)*info.nvar+1:info.nlag*info.nvar,:)=[A{info.nlag}/A0 info.extraF];



            activity_index = zeros(12,1);
            oil_price      = zeros(12,1);

            for h=0:1:11

                tmp = (info.J'*((info.bigF')^h)*info.J)*L0old;

                activity_index(h+1) =tmp(2,1);
                oil_price(h+1)      =tmp(3,1);

             end


           label = all(activity_index<0)*all(oil_price>0);



    end
    


    St_old = label;







