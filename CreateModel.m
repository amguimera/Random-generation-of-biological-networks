function [modelObj] = CreateModel(ModelName, InactiveEnzymes, ActiveEnzymes,InitialAbundancesInactiveEnzymes, InitialAbundancesActiveEnzymes,RateConstantValues, RateConstantNames, Reactions, ReactionRates, AcuteSignalStrength, Timecourse_Duration, EdgeTabulation)

%Create model from the randomly generated reactions

%Define model object
modelObj = sbiomodel (ModelName); %Define model name

%Create compartment
compObj = addcompartment(modelObj, 'cell');

%Create all inactive enzyme species
for i=1:numel(InactiveEnzymes)
    species_name=cell2mat(InactiveEnzymes(i));
    addspecies (compObj, species_name,InitialAbundancesInactiveEnzymes);
end

%Create all active enzyme species
for i=1:numel(ActiveEnzymes)
    species_name=cell2mat(ActiveEnzymes(i));
    addspecies (compObj, species_name,InitialAbundancesActiveEnzymes);
end

%Set rate constants
for i=1:numel(RateConstantValues)
    constant_name=cell2mat(RateConstantNames(i));
    constant_val=RateConstantValues(i);
    addparameter(modelObj, constant_name, constant_val,'ConstantValue',false);
end

% Add the reactions
for i=1:numel(Reactions)
    reaction=cell2mat(Reactions(i));
    rate=cell2mat(ReactionRates(i));
    addreaction(modelObj, reaction, 'ReactionRate', rate);
end

%Delete unused model components

unused = findUnusedComponents(modelObj);
if numel(unused)>0
    Unused={};
    for i=1:numel(unused)
        Unused{end+1}=unused(i).Name;
    end
    for i=numel(modelObj.Species):-1:1
        for ii=1:numel(Unused)
            if strcmp(modelObj.Species(i).Name,str2mat(Unused(ii)))
                delete(modelObj.Species(i))
                break  %since there should only be one occurrence
            end
        end
    end
end


%Create an input reaction that will model a stimulus

%Choose scale-free random species to be affected by input and a rate constant
%{
Accept_Node=0;
while Accept_Node==0
    RandomNode=rand(1); %Choose random number
    Rndm_Species=RandomSelection(EdgeTabulation,RandomNode); %Check what node it corresponds to (Scale free)
    Random_Species=['E' num2str(Rndm_Species)];
    for i=1:numel(modelObj.Species)
        species=modelObj.Species(i).Name;
        if strcmp(species,Random_Species)==1
            Accept_Node=1;
            break
        end
    end
end
%}


%Choose a uniformly random species to be affected by input and a rate constant
UsedSpecies=1:numel(modelObj.Species);
Random_Species='Fx'; 
while Random_Species(1)~='E' %Species affected by input can only be an inactive enzyme
    Rndm_Species=datasample(UsedSpecies,1);
    Random_Species=modelObj.Species(Rndm_Species).Name;
end
 

Random_Species_F_equivalent=['F' Random_Species(2:end)];
%Rndm_RateConstant=datasample(RateConstantValues,1);
Fastest_RateConstant=max(RateConstantValues); %Use fastest rate constant
constant_name='InputConstant';
input_name='SignallingInput';
reaction=[input_name ' + ' Random_Species  ' -> ' Random_Species_F_equivalent];
rate=[input_name '*' constant_name '*' Random_Species];
addspecies (compObj, input_name,0);
addparameter(modelObj, constant_name, Fastest_RateConstant,'ConstantValue',false);
addreaction(modelObj, reaction, 'ReactionRate', rate);

%Add event
EventTime=round(Timecourse_Duration*0.5);
Timing=['time>=' num2str(EventTime)];
Input=[input_name '=' num2str(AcuteSignalStrength)];
event1 = addevent(modelObj,Timing,Input);


end
