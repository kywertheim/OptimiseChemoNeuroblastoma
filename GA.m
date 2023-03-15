% This code implements a genetic algorithm (GA) to find 12 optimal chemotherapy schedules
% with one to 12 cycles, given a neuroblastoma with a specific initial clonal composition.
% The 12 solutions can be improved upon using the fmincon function, which
% searches the vicinity of each solution proposed by the GA.

%%%%%%%%%%%% Step 0: Empty the workspace. %%%%%%%%%%%%
clear all

%%%%%%%%%%%% Step 1: Parameterise clonal evolution. %%%%%%%%%%%%
c=1;    % Importance of interclonal competition, relative to intraclonal competition. 1 means they are equally important.
K=1e10; % Carrying capacity of the tumour.
mu=1e-4;  % Rate of acquiring a random mutation, thus allowing interclonal competition to drive evolution.


%%%%%%%%%%%% Step 2: Parameterise clonal growth. %%%%%%%%%%%%
r_00=0.008506304804128;     % Maximum growth rate of the fully sensitive clone.

% This block defines the relationship between the fully sensitive, the
% mildly resistant, and the strongly resistant clones.
alpha_00=0.031733277051719;    
alpha_01=0.028724723932730;                   
alpha_02=0.027908502236813;                   
a=alpha_01/alpha_00;        % Scaling factor, fully sensitive to mildy resistant.
b=alpha_02/alpha_01;        % Scaling factor, mildy resistant to strongly resistant.
r_01=r_00*a;
r_02=r_01*b;

% This finalises the nine maximum clonal growth rates.
r1=r_00;
r2=r_01;
r3=r_02;
r4=r_01;
r5=r4*a;
r6=r5*b;
r7=r_02;
r8=r7*a;
r9=r8*b;
             
%%%%%%%%%%%% Step 3: Parameterise phenotypic adaptation. %%%%%%%%%%%%
phi2_max=2;       % Maximum level of phenotypic adaptation.
tau_max=10;       % Maximum memory period associated with phenotypic adaptation. A clone will neglect a dose that was present beyond this period.
drug_c_th=0.01;   % Minimum drug concentration needed to trigger phenotypic adaptation.

%%%%%%%%%%%% Step 4: Parameterise vincristine's cytotoxic activity. %%%%%%%%%%%%
alpha_1=1.868604957378585e+05;
beta_1=0.6690469279;
m1_00=2681.9104898293; % Maximum cytotoxic effect on sensitive cells [cell/(ng/mL)].
m1_01=450; % Maximum cytotoxic effect on mildly resistant (VCR-10) cells [cell/(ng/mL)].
m1_02=400; % Maximum cytotoxic effect on strongly resistant (VCR-20) cells [cell/(ng/mL)].

%%%%%%%%%%%% Step 5: Parameterise cyclophosphamide's cytotoxic activity. %%%%%%%%%%%%
alpha_2 = 0.007458099157534;
m2_00=0.000795477343442; % Maximum cytotoxic effect on sensitive cells [cell/(micro M)].
m2_01=0.5*m2_00; % Maximum cytotoxic effect on mildly resistant (CPM-16) cells [cell/(micro M)].
m2_02=0.3*m2_00; % Maximum cytotoxic effect on strongly resistant (CPM-32) cells [cell/(micro M)].
CPM_molar_mass=261.086; % Molar mass of cyclophosphamide [g/mol].

%%%%%%%%%%%% Step 6: Parameterise the patient's individual properties. %%%%%%%%%%%%
% This block describes a three-year-old child whose height is 80 cm and weight is 15 kg.
h_child=80;
kg_child=15;
blood_vol_child=1.125;  % Volume of blood in the child [L] (https://linkprotect.cudasvc.com/url?a=https%3a%2f%2fwww.omnicalculator.com%2fhealth%2fpediatric-blood-volume&c=E,1,wPwTWuzMGzRsLJwrbeUfd6jZdSlTx1FBnX_9sH-BxCG-or0NVEppUQKl7B2hPNakMWHpxCHZUMjfNlU8u4SB5R_wYxv9i5c5R5cZo1tpxdVdsHNGEKLS4vw9j31L&typo=1).
TBW=0.58*kg_child;      % Volume of water in the child.
BSA = sqrt(h_child*kg_child/3600);  % Body surface area of the child [m2].

% This block calculates the clearance rates with respect to the two drugs in the child.
CL_VCR_child=0.228*60; % VCR clearance rate [L/h/m2].
CL_CPM_child=1.77;   % CPM clearance rate [L/h/m2].
z_CPM_child=CL_CPM_child*BSA/TBW;
z_VCR_child=CL_VCR_child*BSA/TBW;

%%%%%%%%%%%% Step 7: Define the chemotherapy schedule. %%%%%%%%%%%%
max_n_cycles=12;    % Maximum number of cycles.
length_cycle=14;    % Duration of one cycle in days.
max_dosage=2; % Maximum tolerated dose. Same value but different units for both VCR and CPM.

%%%%%%%%%%%% Step 8: Configure the genetic algorithm. %%%%%%%%%%%%
N_iterations = 1000; % Number of generations.
N_chromosomes = 100; % Number of chromosomes in each generation.
Probability_mutation = 0.05; % Probability of mutation.
best_chromosome = zeros(N_iterations, 2*max_n_cycles, max_n_cycles); % 3D matrix to save the best chromosomes. If the number of cycles is smaller than max_n_cycles, the last 2*(max_n_cycles-N) genes in a chromosome will take the value of zero.
max_fitness = zeros(N_iterations, max_n_cycles); % 2D matrix to save the fitness of the best chromosome in each generation and for each number of cycles.

%%%%%%%%%%%% Step 9: Apply the genetic algorithm to a broad range of initial compositions. %%%%%%%%%%%%
for res=5:5:25 % Initial percentage of the cancer cells that are resistant.
    for heterogeneity=1:7 % This variable indicates the resistant clones that are present at the start of the simulation.
    % 1: mildly VCR-resistant clone only.
    % 2: mildly VCR-resistant clone and mildly CPM-resistant clone only.
    % 3: mildly CPM-resistant clone only.
    % 4: mildly VCR-resistant clone, mildly CPM-resistant clone, strongly VCR-resistant clone, and strongly CPM-resistant clone only.
    % 5: strongly VCR-resistant clone only.
    % 6: strongly CPM-resistant clone only.
    % 7: strongly VCR-resistant clone and strongly CPM-resistant clone only.
        for cycle=1:max_n_cycles  % Experiment with the number of chemotherapy cycles.
            control_var=2*cycle; % Number of control variables (VCR and CPM doses in each cyle)
            %%%%%%%%%%%% Step 9a: Create the initial generation of chromosomes. %%%%%%%%%%%%
            chromosomes = max_dosage*rand(N_chromosomes, control_var); %Each chromosome has control_var genes. The first half are VCR doses and the other half are CPM doses.
            %Let the chromosomes evolve.
            for iteration=1:N_iterations
                disp(res)
                disp(heterogeneity)
                disp(cycle)
                disp(iteration)
                %%%%%%%%%%%% Step 9b: Evaluate each chromosome's fitness in the current generation. %%%%%%%%%%%%
                fitness = zeros(N_chromosomes, 1);
                for chromosome = 1:N_chromosomes
                    VCR_dosage=chromosomes(chromosome,1: control_var/2)*1e6/(TBW*1e3)*BSA;% Control variables, VCR doses in the right physical units [ng/mL].
                    CPM_dosage=chromosomes(chromosome, control_var/2+1:end)*BSA/CPM_molar_mass/TBW*1e6; % Control variables, CPM doses in the right physical units [microM].
                    
                    % Combine the parameters.
                    pars(1)=K;
                    pars(2)=c;
                    pars(3)=r1;
                    pars(4)=r2;
                    pars(5)=r3;
                    pars(6)=r4;
                    pars(7)=r5;
                    pars(8)=r6;
                    pars(9)=r7;
                    pars(10)=r8;
                    pars(11)=r9;
                    pars(12)=mu;
                    pars(13)=m1_00;
                    pars(14)=m1_01;
                    pars(15)=m1_02;
                    pars(16)=alpha_1;
                    pars(17)=beta_1;
                    pars(18)=m2_00;
                    pars(19)=m2_01;
                    pars(20)=m2_02;
                    pars(21)=alpha_2;
                    pars(22)=phi2_max;
                    pars(23)=tau_max;
                    pars(24)=drug_c_th;
                    pars(25)=z_CPM_child;
                    pars(29)=z_VCR_child;
                    pars(26)=VCR_dosage(1); % VCR dose in the first chemotherapeutic cycle.
                    pars(27)=CPM_dosage(1); % CPM dose in the first chemotherapeutic cycle.
                    pars(28)=length_cycle*24; % Turn days to hours.

                    % Initialise the tumour.                   
                    N1_0=K/2;% Initial total cell count.
                    perc_res=res; % Initial percentage of the cancer cells that are resistant.                    
                    X_iniz=[N1_0*(100-perc_res)/100 0 0 0 0 0 0 0 0 0 0 0];% Initial number of fully sensitive cells.
                    clones=heterogeneity; % This variable indicates the resistant clones that are present at the start of the simulation.
                    % This block distributes the resistant cells between the resistant clones
                    % present at the start of the simulation.                    
                    if clones==1
                        X_iniz(2)=N1_0*perc_res/100; % mildly VCR-resistant clone.
                    elseif clones==2
                        X_iniz(2)=N1_0*perc_res/100/2; % mildly VCR-resistant clone.
                        X_iniz(4)=N1_0*perc_res/100/2; % mildly CPM-resistant clone.
                    elseif clones==3
                        X_iniz(4)=N1_0*perc_res/100; % mildly CPM-resistant clone.
                    elseif clones==4
                        X_iniz(2)=N1_0*perc_res/100*1/3; % mildly VCR-resistant clone.
                        X_iniz(4)=N1_0*perc_res/100*1/3; % mildly CPM-resistant clone.
                        X_iniz(3)=N1_0*perc_res/100*1/6; % strongly VCR-resistant clone.
                        X_iniz(7)=N1_0*perc_res/100*1/6; % strongly CPM-resistant clone.
                     elseif clones==5
                        X_iniz(3)=N1_0*perc_res/100; % strongly VCR-resistant clone.
                     elseif clones==6
                        X_iniz(7)=N1_0*perc_res/100; % strongly CPM-resistant clone.
                     elseif clones==7
                        X_iniz(3)=N1_0*perc_res/100/2; % strongly VCR-resistant clone.
                        X_iniz(7)=N1_0*perc_res/100/2; % strongly CPM-resistant clone.
                    end

                    x0=X_iniz; % Summarise the initial cell counts.

                    % Set up the ordinary differential equations governing the system.
                    odefun = @(t,x) logistic_growth(t,x,pars);
                    options1 = odeset('Refine',1);
                    options = odeset(options1,'NonNegative',1:12);
                    
                    % This block simulates the first cycle of the chemotherapy schedule.
                    t_exp=length_cycle*24; % Specify the extent in the temporal domain: one cycle in hours.
                    [T,X] = ode45(odefun,[0,t_exp],x0,options); % Use a numerical solver to simulate the population dynamics in the first chemotherapeutic cycle.
                    
                    % This block simulates the remaining cycles one by one.
                    for i=1:cycle-1
                        x0=X(end,:);   % The final conditions of the previous cycle are the initial conditions of the current cycle.
                        pars(26)=VCR_dosage(i+1); % VCR dose in the current cycle.
                        pars(27)=CPM_dosage(i+1); % CPM dose in the current cycle.
                        odefun = @(t,x) logistic_growth(t,x,pars);
                        [T_new,X_new] = ode45(odefun,[t_exp*i,t_exp*(i+1)],x0,options);
                        T=[T;T_new]; %Add the new results to the overall results.
                        X=[X;X_new]; %Add the new results to the overall results.
                    end

                    %The chromosome's fitness is the difference between the
                    %carrying capacity and the final population size.
                    fitness(chromosome,1)=K-sum(X(end,1:9));        
                end
                % Record the best chromosome and its fitness in this generation.
                [dummy_maxfit, dummy_maxfitpos] = max(fitness);
                max_fitness(iteration,cycle) = dummy_maxfit;
                best_chromosome(iteration,1:control_var,cycle) = chromosomes(dummy_maxfitpos,:);
                
                % Check that the maximum fitness in this generation is at least equal to that in the last generation.
                disp(dummy_maxfit)
                if iteration > 1
                    if dummy_maxfit < max_fitness(iteration-1,cycle)
                       disp('Maximum fitness has decreased!')
                       break
                    end
                end

                %%%%%%%%%%%% Step 9c: Create the next generation of chromosomes. %%%%%%%%%%%%
                % Select the top half of the current generation and calculate their mating probabilities.
                [chromosomes_offspring, probabilities_mating] = selection(chromosomes, fitness, N_chromosomes,control_var);

                % Pair up the selected chromosomes and let them cross over to create half of the next generation.
                % The top 10 chromosomes in the current generation are allowed to enter the next generation without further consideration.
                % Create 40 random chromosomes to complete the next generation.
                new_chromosomes=ceil(N_chromosomes*4/10); % Specify the number of random chromosomes to be added.
                chromosomes_offspring = crossover(chromosomes_offspring, N_chromosomes, probabilities_mating,max_dosage,control_var,new_chromosomes);

                %%%%%%%%%%%% Step 9d: Mutation. %%%%%%%%%%%%
                % Attempt to mutate the 50 daughter chromosomes only.
                chromosomes_offspring = mutation(chromosomes_offspring, N_chromosomes, Probability_mutation, max_dosage,new_chromosomes,control_var);

                %%%%%%%%%%%% Step 9e: Enter the next generation. %%%%%%%%%%%%
                % Confirm the new population of chromosomes for the next generation.
                chromosomes = chromosomes_offspring;
            end
        end
        % Save the results for this initial composition.
        filename = join(['GA_results_res', num2str(perc_res),'_clonecomp', num2str(clones)]);
        save(filename)
    end
end

%%%%%%%%%%%% Ordinary differential equations governing the system. %%%%%%%%%%%%
function dNdt = logistic_growth(t,x,pars)
    
    % This block parameterises the equations.
    K=pars(1);
    c=pars(2);
    r00=pars(3);
    r01=pars(4);
    r02=pars(5);
    r10=pars(6);
    r11=pars(7);
    r12=pars(8);
    r20=pars(9);
    r21=pars(10);
    r22=pars(11);
    mu=pars(12);
    d_v_0=pars(13);
    d_v_1=pars(14);
    d_v_2=pars(15);
    delta_1=pars(16);
    gamma_1=pars(17);
    d_c_0=pars(18);
    d_c_1=pars(19);
    d_c_2=pars(20);
    delta_2=pars(21);
    R_max=pars(22);
    gg_pr=pars(23);
    drug_c_th=pars(24);
    z_CPM=pars(25);
    z_VCR=pars(29);

    % This block defines the two dosages in the current cycle.
    dose_VCR=pars(26); % Amount of VCR.
    dose_CPM=pars(27); % Amount of CPM.
    freq_dose=pars(28); % Cycle duration.
    h_admin_VCR=2*24; % This indicates when in the cycle is VCR administered.
    h_admin_CPM=4*24; % This indicates when in the cycle is CPM administered.    

    % This block initialises the variables.
    N00=x(1);
    N01=x(2);
    N02=x(3);
    N10=x(4);
    N11=x(5);
    N12=x(6);
    N20=x(7);
    N21=x(8);
    N22=x(9);
    y1=x(10); % VCR concentration.
    y2=x(11); % CPM concentration.    
    tau=x(12);% Extent to which phenotypic adaptation is active.
    
    % This block simplifies the algebra.
    reduced_proliferation=(1+R_max*(tau/gg_pr)/R_max);
    reduce_drug_c=(1+R_max*(tau/gg_pr));

    % This block defines the ordinary differential equations.
    dNdt=[
    %dot N00    
    r00/reduced_proliferation*(1-(N00+c*(N01+N02+N10+N11+N12+N20+N21+N22))/K)*N00*(1-2*mu)+mu*(r01/reduced_proliferation*(1-(N01+c*(N00+N02+N10+N11+N12+N20+N21+N22))/K)*N01+r10/reduced_proliferation*(1-(N10+c*(N00+N01+N02+N11+N12+N20+N21+N22))/K)*N10)-(y1>0)*d_v_0/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N00-d_c_0/(1+delta_2*y2)*y2/reduce_drug_c*N00;
    %dot N01    
    r01/reduced_proliferation*(1-(N01+c*(N00+N02+N10+N11+N12+N20+N21+N22))/K)*N01*(1-3*mu)+mu*(r00/reduced_proliferation*(1-(N00+c*(N01+N02+N10+N11+N12+N20+N21+N22))/K)*N00+r02/reduced_proliferation*(1-(N02+c*(N00+N01+N10+N11+N12+N20+N21+N22))/K)*N02+r11/reduced_proliferation*(1-(N11+c*(N00+N01+N02+N10+N12+N20+N21+N22))/K)*N11)-(y1>0)*d_v_1/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N01-d_c_0/(1+delta_2*y2)*y2/reduce_drug_c*N01;
    %dot N02    
    r02/reduced_proliferation*(1-(N02+c*(N00+N01+N10+N11+N12+N20+N21+N22))/K)*N02*(1-2*mu)+mu*(r01/reduced_proliferation*(1-(N01+c*(N00+N02+N10+N11+N12+N20+N21+N22))/K)*N01+r12/reduced_proliferation*(1-(N12+c*(N00+N01+N02+N10+N11+N20+N21+N22))/K)*N12)-(y1>0)*d_v_2/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N02-d_c_0/(1+delta_2*y2)*y2/reduce_drug_c*N02;
    %dot N10    
    r10/reduced_proliferation*(1-(N10+c*(N00+N01+N02+N11+N12+N20+N21+N22))/K)*N10*(1-3*mu)+mu*(r00/reduced_proliferation*(1-(N00+c*(N01+N02+N10+N11+N12+N20+N21+N22))/K)*N00+r20/reduced_proliferation*(1-(N20+c*(N00+N01+N02+N10+N11+N12+N21+N22))/K)*N20+r11/reduced_proliferation*(1-(N11+c*(N00+N01+N02+N10+N12+N20+N21+N22))/K)*N11)-(y1>0)*d_v_0/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N10-d_c_1/(1+delta_2*y2)*y2/reduce_drug_c*N10;
    %dot N11    
    r11/reduced_proliferation*(1-(N11+c*(N00+N01+N02+N10+N12+N20+N21+N22))/K)*N11*(1-4*mu)+mu*(r01/reduced_proliferation*(1-(N01+c*(N00+N02+N10+N11+N12+N20+N21+N22))/K)*N01+r10/reduced_proliferation*(1-(N10+c*(N00+N01+N02+N11+N12+N20+N21+N22))/K)*N10+r12/reduced_proliferation*(1-(N12+c*(N00+N01+N02+N10+N11+N20+N21+N22))/K)*N12+r21/reduced_proliferation*(1-(N21+c*(N00+N01+N02+N10+N11+N12+N20+N22))/K)*N21)-(y1>0)*d_v_1/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N11-d_c_1/(1+delta_2*y2)*y2/reduce_drug_c*N11;
    %dot N12    
    r12/reduced_proliferation*(1-(N12+c*(N00+N01+N02+N10+N11+N20+N21+N22))/K)*N12*(1-3*mu)+mu*(r11/reduced_proliferation*(1-(N11+c*(N00+N01+N02+N10+N12+N20+N21+N22))/K)*N11+r02/reduced_proliferation*(1-(N02+c*(N00+N01+N12+N10+N11+N20+N21+N22))/K)*N12+r22/reduced_proliferation*(1-(N22+c*(N00+N01+N02+N10+N11+N12+N20+N21))/K)*N22)-(y1>0)*d_v_2/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N12-d_c_1/(1+delta_2*y2)*y2/reduce_drug_c*N12;
    %dot N20    
    r20/reduced_proliferation*(1-(N20+c*(N00+N01+N02+N10+N11+N12+N21+N22))/K)*N20*(1-2*mu)+mu*(r10/reduced_proliferation*(1-(N10+c*(N00+N01+N02+N11+N12+N20+N21+N22))/K)*N10+r21/reduced_proliferation*(1-(N21+c*(N00+N01+N02+N10+N11+N12+N20+N22))/K)*N21)-(y1>0)*d_v_0/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N20-d_c_2/(1+delta_2*y2)*y2/reduce_drug_c*N20;
    %dot N21    
    r21/reduced_proliferation*(1-(N21+c*(N00+N01+N02+N10+N11+N12+N20+N22))/K)*N21*(1-3*mu)+mu*(r20/reduced_proliferation*(1-(N20+c*(N00+N01+N02+N10+N11+N12+N21+N22))/K)*N20+r11/reduced_proliferation*(1-(N11+c*(N00+N01+N02+N10+N12+N20+N21+N22))/K)*N11+r22/reduced_proliferation*(1-(N22+c*(N00+N01+N02+N10+N11+N12+N20+N21))/K)*N22)-(y1>0)*d_v_1/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N21-d_c_2/(1+delta_2*y2)*y2/reduce_drug_c*N21;
    %dot N22    
    r22/reduced_proliferation*(1-(N22+c*(N00+N01+N02+N10+N11+N12+N20+N21))/K)*N22*(1-2*mu)+mu*(r21/reduced_proliferation*(1-(N21+c*(N00+N01+N02+N10+N11+N12+N20+N22))/K)*N21+r12/reduced_proliferation*(1-(N12+c*(N00+N01+N02+N10+N11+N20+N21+N22))/K)*N12)-(y1>0)*d_v_2/(1+delta_1*y1^gamma_1)*y1/reduce_drug_c*N22-d_c_2/(1+delta_2*y2)*y2/reduce_drug_c*N22;
    %dot y1 --> VCR drug concentration
    dose_VCR/h_admin_VCR*(rem(t,freq_dose)<h_admin_VCR)-z_VCR*y1; 
    %dot y2 --> CPM drug concentration   
    dose_CPM/h_admin_CPM*(rem(t,freq_dose)<h_admin_CPM)-z_CPM*y2; 
    % days under drug stress
    1/24*(y1+y2>drug_c_th)*(tau<gg_pr)-1/24*(y1+y2<drug_c_th)*(tau>0);
    ];
end

%%%%%%%%%%%% Selection of the fittest chromosomes in a generation. %%%%%%%%%%%%
function [chromosomes_offspring, probabilities_mating] = selection(chromosomes, fitness, N_chromosomes,control_var)
%Identify the chromosomes in the top half of the generation,
%add them to the next generation, calculate their mating probabilities, and
%turn the probabilities into cumulative probabilities.

    %Find out the top half of the chromosomes in terms of their fitness values and indices.
    [dummy_maxfit, dummy_maxfitpos] = maxk(fitness, N_chromosomes/2);
    
    % Create data structures for recording.
    chromosomes_offspring = zeros(N_chromosomes, control_var);
    probabilities_mating = zeros(N_chromosomes/2, 1);

    % Extract the chromosomes in the top half and calculate their mating probabilities.
    for i=1:N_chromosomes/2
        probabilities_mating(i)=dummy_maxfit(i)/sum(dummy_maxfit); % Turn the fitness value of the current chromosome into a mating probability.
        j=dummy_maxfitpos(i); % Locate the current chromosome.
        chromosomes_offspring(i, :) = chromosomes(j, :); % Add the selected chromosome to the next generation.
    end

    % Turn the mating probabilities into cumulative probabilities.
    for i=2:N_chromosomes/2
        probabilities_mating(i)=probabilities_mating(i-1)+probabilities_mating(i);
    end

end

%%%%%%%%%%%% Creation of the next generation of chromosomes. %%%%%%%%%%%%
function chromosomes_offspring_copy = crossover(chromosomes_offspring, N_chromosomes, probabilities_mating, max_VCR_intermittent,control_var,new_chromosomes)
chromosomes_offspring_copy=chromosomes_offspring;
%After selection, let the top 10 chromosomes enter the next generation
%without further consideration, let the 50 selected chromosomes reproduce
%50 daughter chromosomes, and randomly generate 40 new chromosomes
%(migration).

    for i=1:N_chromosomes/4
    % Let the top half of the current generation exchange genes to create
    % 50 daughter chromosomes for the next generation.
    
       % Pick two chromosomes from the top half of the current generation
       % according to the mating probabilities (cumulative form).
       dummy1 = rand;
       dummy2 = rand;
       index1 = 0;
       index2 = 0;
       for j=1:N_chromosomes/2
          if dummy1<probabilities_mating(j)
              index1 = j;
              break
          end
       end
       for j=1:N_chromosomes/2
          if dummy2<probabilities_mating(j)
              index2 = j;
              break
          end
       end
       
       % Let the two mating chromosomes exchange genes to produce daughter chromosomes for the next generation.
       index_crossover = randi([1 control_var]); % Crossover point.
       beta = rand; % Extent of crossover.
       for j=1:control_var
           if j<=index_crossover % If the gene is at or before the crossover point, let the chromosomes exchange genes before reproducing.
               chromosomes_offspring_copy(N_chromosomes/2+2*i-1-new_chromosomes, j) = (1-beta)*chromosomes_offspring(index1, j)+beta*chromosomes_offspring(index2, j);
               chromosomes_offspring_copy(N_chromosomes/2+2*i-new_chromosomes, j) = (1-beta)*chromosomes_offspring(index2, j)+beta*chromosomes_offspring(index1, j);
           else % Otherwise, the chromosomes will just pass on this pair of genes directly without exchanging.
               chromosomes_offspring_copy(N_chromosomes/2+2*i-1-new_chromosomes, j) = chromosomes_offspring(index1, j);
               chromosomes_offspring_copy(N_chromosomes/2+2*i-new_chromosomes, j) = chromosomes_offspring(index2, j);        
           end
       end
    end

    % Create random chromosomes.
    chromosomes_offspring_copy(end-new_chromosomes+1:end,:)=max_VCR_intermittent*rand(new_chromosomes, control_var);
end

%%%%%%%%%%%% Mutate the daughter chromosomes produced by the selected chromosomes. %%%%%%%%%%%%
function chromosomes_offspring = mutation(chromosomes_offspring, N_chromosomes, Probability_mutation, max_VCR_intermittent,new_chromosomes,control_var)
% Try to mutate the daughter chromosomes, not the 10 best chromosomes or
% the 40 randomly created chromosomes.
    for i=N_chromosomes/2-new_chromosomes+1:N_chromosomes-new_chromosomes
        for j=1:control_var
           if rand < Probability_mutation
               chromosomes_offspring(i, j) = max_VCR_intermittent*rand; % Mutation means picking a random non-zero number smaller than the maximum dose.
           end
        end
    end
end