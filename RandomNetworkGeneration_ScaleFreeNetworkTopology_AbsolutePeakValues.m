%%Random network generation to test robustness of biological systems to
%%constitutive signals

%%Developed by Alvaro Martinez Guimera July 2019


%Shuffle random number generator seeding
rng('shuffle')

%User-inputs
NumOfNetworks=10000; %Number of networks to be generated
NumOfNodes=25;  %Define number of nodes (enzymes) - note that this number will be duplicated
RateConstantScale=0.01; %defines scale of speed of enzymatic reactions. This number will be converted into a distribution to be sampled from.
InitialAbundanceInactiveEnzymes=10;
InitialAbundanceActiveEnzymes=10;
AcuteSignalStrength=InitialAbundanceInactiveEnzymes*10;  %Magnitude of network stimulus
FractionOfPossibleEdges=0.04; %Fraction of interaction combinations explored

%Generate Enzyme identifiers and derive number on interactions to include in the network
EnzymeIdentifiers=1:NumOfNodes;
InteractionCombinations=((NumOfNodes^2)-NumOfNodes)*2; %Squaring due to pairwise interactions
%between nodes. The subtraction of the number of nodes is due to
%self-regulation not being allowed. Multiplied by two because there are two
%sets of pairwise interactions (F -> E) and (F -> F);
Interactions=round(FractionOfPossibleEdges*InteractionCombinations);

%Since rates will be modelled as mass action, the
%number of constants should match the number of interactions
RateConstantNames={};
for i=1:Interactions
    RateConstantNames{end+1}=['k' num2str(i)];
end

%Define population of inactive enzymes (E) and active enzymes (F). 
E={}; %Identifiers for inactive enzymes
F={}; %Identifies for active enzymes
Ei=InitialAbundanceInactiveEnzymes; %Initial abundance of inactive enzymes
Fi=InitialAbundanceActiveEnzymes; %Initial abundance of active enzymes
for i=1:NumOfNodes
    E{end+1}=['E' num2str(i)];
    F{end+1}=['F' num2str(i)];
end

%Keep track of models that satisfy requisites
Models_saved=1;

for n=1:NumOfNetworks
    
    %Reaction string storage
    Reactions={};
    ReactionRates={};
    
    %Rate constant values
    RateConstantValues=[]; %Get appended as reactions are generated 

    %Create tabulation table to bias probability of node selection for a scale-free network
    Counters=(zeros(NumOfNodes,1)+1);
    EdgeTabulation=[EnzymeIdentifiers' Counters];
    NormalisationFactor=1/sum(EdgeTabulation(:,2)); %Normalise tabulation table probabilities
    EdgeTabulation(:,2)=EdgeTabulation(:,2)*NormalisationFactor;

    
    for i=1:Interactions %For a fraction of the total number of combinations

        %Choose an initial random 'F' 
        RandomNode=rand(1);
        Rndm_F = RandomSelection(EdgeTabulation,RandomNode);
        
        %Update the tabulation table
        Counters(Rndm_F)=Counters(Rndm_F)+1; %Update tabulation table
        NormalisationFactor=1/sum(Counters); %Renormalise tabulation table probabilities
        Updated_Probabilities=Counters*NormalisationFactor;
        EdgeTabulation(:,2)=Updated_Probabilities;

        %Choose a random binary outcome for another random node  
        EdgeType=rand(1);
        RandomNode2=rand(1);
        if EdgeType<=0.5 % 'F' acts on another 'F' (deactivation)
            Rndm_F_2=RandomSelection(EdgeTabulation,RandomNode2);
            while Rndm_F_2 == Rndm_F %No self-inhibition allowed
                RandomNode2=rand(1);
                Rndm_F_2=RandomSelection(EdgeTabulation,RandomNode2); %Select another node
            end
            React=[cell2mat(F(Rndm_F)) ' + ' cell2mat(F(Rndm_F_2)) ' -> ' cell2mat(E(Rndm_F_2)) ' + ' cell2mat(F(Rndm_F))];
            ReactRate=[cell2mat(F(Rndm_F)) '*' cell2mat(F(Rndm_F_2)) '*' cell2mat(RateConstantNames(i))];
            Reactions{end+1}=React; %Update reaction list
            ReactionRates{end+1}=ReactRate; %Update reaction rate list
            RateConstantValues=[RateConstantValues; RateConstantScale*2]; %Deactivation rates are assumed to be double that of activation rates
        else
            Rndm_E=RandomSelection(EdgeTabulation,RandomNode2); % 'F' acts on an 'E' (activation)
            while Rndm_E == Rndm_F %No self-activation allowed
                RandomNode2=rand(1);
                Rndm_E=RandomSelection(EdgeTabulation,RandomNode2); %Select another node
            end
            React=[cell2mat(F(Rndm_F)) ' + ' cell2mat(E(Rndm_E)) ' -> ' cell2mat(F(Rndm_E)) ' + ' cell2mat(F(Rndm_F))];
            ReactRate=[cell2mat(F(Rndm_F)) '*' cell2mat(E(Rndm_E)) '*' cell2mat(RateConstantNames(i))];
            Reactions{end+1}=React; %Update reaction list
            ReactionRates{end+1}=ReactRate; %Update reaction rate list
            RateConstantValues=[RateConstantValues; RateConstantScale]; %Deactivation rates are assumed to be double that of activation rates
        end
    end

    %Add variation to rate constants and cap those that become negative to 1% of origina parameter value 
    mean=1; %Average value for rate constant multiplier is 1
    variance=mean/10;
    variation = sqrt(variance)*randn(Interactions,1)+mean;
    idxs=find(variation<=0);
    if numel(idxs)>0
        for i=1:numel(idxs)
            idx=idxs(i);
            variation(idx)=RateConstantScale/100; % 1% of original value is the minimum rate constant value 
        end
    end
    RateConstantValues = RateConstantValues.*variation;
    
    %Define time course duration
    Timecourse_Duration=(1/RateConstantScale)*2; %One for equilibration and one for simulation
    
    %Create Simbiology model
    ModelName=['RandomNetwork_' num2str(n)];
    [modelObj] = CreateModel(ModelName, E, F, Ei, Fi, RateConstantValues, RateConstantNames, Reactions, ReactionRates, AcuteSignalStrength, Timecourse_Duration, EdgeTabulation);
    
    %Perform deterministic simulation
    cs = getconfigset(modelObj,'active');
    cs.SolverType = 'ode15s';  %ODE solver - ODE45 non-stiff, ODE15s & ODE23 = stiff
    cs.SolverOptions.AbsoluteTolerance= 1.0e-12;
    cs.SolverOptions.RelativeTolerance= 1.0e-6;
    cs.StopTime = Timecourse_Duration; %Simulation stop time
    cs.CompileOptions.UnitConversion = false;  %No unit conversion
    Simulation_Output=sbiosimulate(modelObj); %Simulate model
    %sbioplot(Simulation_Output);  %plot simulation output
    
    %Define relevant data range in species and time
    StimulatedSpeciesThreshold=10;
    EventTime=round(Timecourse_Duration/2);
    [number, MinTimeIndex ] = min( abs( Simulation_Output.Time-EventTime ) );%Identify index closes to Event time in order not to include peaks/troughs arising from system equilibration
    F_index = find(contains(Simulation_Output.DataNames,'F')); %Find index range of F species
     
    %Look for threshold number of stimulated F species
    StimulatedSpecies=0; %counter
    Peak_data=[];
    for i=F_index(1):F_index(end) %for each simulated F species
        
        %Check for peaks
        RelevantData=Simulation_Output.Data(MinTimeIndex:end,i);
        [PeakMagnitude,location,width,prominence]=findpeaks(RelevantData,'MinPeakProminence',1); %Minimum peak magnitude of 1 a.u
        
        %Check for troughs
        RelevantData2=RelevantData*-1; 
        [TroughMagnitude,location2,width2,prominence2]=findpeaks(RelevantData2,'MinPeakProminence',1); %Minimum trough magnitude of 1 a.u
       
        %If both peaks and troughs are found select most prominent
        if numel(PeakMagnitude)>0 && numel(TroughMagnitude)>0 
            Peak=max(prominence);
            Trough=max(prominence2);
            if Peak>Trough
                peak_info=[i max(PeakMagnitude)];
                Peak_data=[Peak_data; peak_info]; %store species index and peak magnitude
                StimulatedSpecies=StimulatedSpecies+1;
            else
                peak_info=[i max(TroughMagnitude)*-1];
                Peak_data=[Peak_data; peak_info]; %store species index and peak magnitude
                StimulatedSpecies=StimulatedSpecies+1;
            end
            
        elseif numel(PeakMagnitude)>0 && numel(TroughMagnitude)==0 %if only peaks found
            Peak=max(PeakMagnitude);
            peak_info=[i Peak];
            Peak_data=[Peak_data; peak_info]; %store species index and peak magnitude
            StimulatedSpecies=StimulatedSpecies+1;
            
        elseif numel(PeakMagnitude)==0 && numel(TroughMagnitude)>0 %if only troughs found
            Trough=max(TroughMagnitude)*-1;
            peak_info=[i Trough];
            Peak_data=[Peak_data; peak_info]; %store species index and peak magnitude
            StimulatedSpecies=StimulatedSpecies+1;
        end
             
    end
    
    %Save model and data if threshold of stimulated species is surpassed
    if StimulatedSpecies>=StimulatedSpeciesThreshold
        ModelName2=['RandomNetwork_' num2str(Models_saved)];
        %Save Model as
        sbiosaveproject(ModelName2, 'modelObj') %Simbiology project
        xlswrite(ModelName2,Peak_data); %Save peak data
        %sbmlexport(modelObj, 'File_Name')  %SBML file
        Models_saved=Models_saved+1;
    end


end


