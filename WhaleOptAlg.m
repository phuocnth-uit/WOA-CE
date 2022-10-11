%_________________________________________________________________________%
%  Whale Optimization Algorithm (WOA) source codes                        %
%_________________________________________________________________________%

%__________________________________________
% fobj              = @YourCostFunction
% dim               = number of your variables
% Max_iteration     = maximum number of generations
% SearchAgents_no   = number of search agents
% lb                = [lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub                = [ub1,ub2,...,ubn] where ubn is the upper bound of variable n


% The Whale Optimization Algorithm
function Leader_pos = WhaleOptAlg(SearchAgents_no,Max_iter,lb,ub,dim,fobj)

% initialize position vector and score for the leader
Leader_pos=zeros(1,dim);
Leader_score=inf;                   %change this to -inf for maximization problems


%Initialize the positions of search agents
Positions = rand(SearchAgents_no, dim).*(ub-lb) + lb;

All_fitness = zeros(1, SearchAgents_no);


% Main loop
for t = 1:Max_iter
    parfor i=1:SearchAgents_no
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=(Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % Calculate objective function for each search agent     
        All_fitness(1,i) = fobj(Positions(i,:));
                
    end
    
    for i=1:SearchAgents_no
        % Update the leader
        fitness= All_fitness(1,i);
        if fitness < Leader_score % Change this to > for maximization problem
            Leader_score = fitness; % Update alpha
            Leader_pos = Positions(i,:);
        end
    end
    
    a=2-t*((2)/Max_iter); % a decreases linearly fron 2 to 0 in Eq. (2.3)
    
    % a2 linearly dicreases from -1 to -2 to calculate t in Eq. (3.12)
    a2=-1+t*((-1)/Max_iter);
    
    % Update the Position of search agents 
    for i=1:SearchAgents_no
        r1=rand(); % r1 is a random number in [0,1]
        r2=rand(); % r2 is a random number in [0,1]
        
        A=2*a*r1-a;  % Eq. (2.3) in the paper
        C=2*r2;      % Eq. (2.4) in the paper
        
        
        b=1;               %  parameters in Eq. (2.5)
        l = 2*rand - 1;   %  parameters in Eq. (2.5)
        
        p = rand();        % p in Eq. (2.6)
        
        for j=1:dim
            
            if p<0.5   
                if abs(A)>=1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand=abs(C*X_rand(j)-Positions(i,j));       % Eq. (2.7)
                    Positions(i,j)=X_rand(j)-A*D_X_rand;            % Eq. (2.8)
                 
                elseif abs(A)<1
                    D_Leader=abs(C*Leader_pos(j)-Positions(i,j));   % Eq. (2.1)
                    Positions(i,j)=Leader_pos(j)-A*D_Leader;        % Eq. (2.2)
                end
                
            elseif p>=0.5
              
                distance2Leader=abs(Leader_pos(j)-Positions(i,j));
                % Eq. (2.5)
                Positions(i,j)=distance2Leader*exp(b.*l).*cos(l.*2*pi)+Leader_pos(j);
                
            end
            
        end
    end

end



