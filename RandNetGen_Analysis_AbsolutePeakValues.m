%%Data analysis of Random network generation to test robustness of biological systems to
%%constitutive signals

%%Developed by Alvaro Martinez Guimera July 2019

%Specify time of model event (stimulus)
EventTime=100;

%Retrieve all files in current directory (where random network generation output should be saved to)
Files=dir('RandomNetwork_*.*');

Results=[];
MetaData=[];
for k=1:2:length(Files) %increases in groups of two since files exist in pairs (one for model and one for peak data)
    
    %Retrieve file names
    DataFileName=Files(k+1).name; %Data files are listed after .sbproj files 
    ModelFileName=Files(k).name;
    
    %Import model and data
    sbioloadproject(ModelFileName);
    ReferencePeakData = xlsread(DataFileName);
    StimulatedSpeciesIndex = ReferencePeakData(:,1); 
    StimulatedSpeciesMagnitude = ReferencePeakData(:,2);
    
    PerturbationData=[]; %store data for peak changes for each parameter
    
    %Perform parameter perturbations
    for p=1:numel(modelObj.parameters)
        
        %Modify parameter
        original_parameter=modelObj.parameters(p).value; %retrieve original parameter
        new_parameter=10*original_parameter; 
        modelObj.parameters(p).value=new_parameter; %update parameter value
        
        %Perform deterministic simulation
        Simulation_Output=sbiosimulate(modelObj); %Simulate model
        %sbioplot(Simulation_Output);  %plot simulation output
        
        %Identify index where simulation reaches event time in order not to include peaks/troughs arising from system equilibration
        [number, MinTimeIndex ] = min( abs( Simulation_Output.Time-EventTime ) );
        
        %Calculate changes in peak stimulation compared to reference simulation
        Average_response_data=[];
        for i=1:numel(StimulatedSpeciesIndex) %for each stimulated F species
            
            %Define relevant data
            RelevantSpeciesIndex=StimulatedSpeciesIndex(i);
            RelevantData=Simulation_Output.Data(MinTimeIndex:end,RelevantSpeciesIndex);
            
            %Check for peaks
            [PeakMagnitude,location,width,prominence]=findpeaks(RelevantData); %without minimum peak magnitude 
            
            %Check for troughs
            RelevantData2=RelevantData*-1; 
            [TroughMagnitude,location2,width2,prominence2]=findpeaks(RelevantData2); %without minimum trough magnitude
       
            %If both peaks and troughs are found select most prominent
            if numel(PeakMagnitude)>0 && numel(TroughMagnitude)>0 
                Peak=max(prominence);
                Trough=max(prominence2);
                if Peak>Trough
                    ResponseChange=max(PeakMagnitude)-StimulatedSpeciesMagnitude(i);
                    Average_response_data=[Average_response_data; ResponseChange];
                else
                    ResponseChange=-((max(TroughMagnitude)*-1)-StimulatedSpeciesMagnitude(i));
                    Average_response_data=[Average_response_data; ResponseChange];                   
                end

            elseif numel(PeakMagnitude)>0 && numel(TroughMagnitude)==0 %if only peaks found
                Peak=max(PeakMagnitude);
                ResponseChange=Peak-StimulatedSpeciesMagnitude(i);
                Average_response_data=[Average_response_data; ResponseChange];
                
            elseif numel(PeakMagnitude)==0 && numel(TroughMagnitude)>0 %if only troughs found
                Trough=max(TroughMagnitude)*-1;
                ResponseChange=-(Trough-StimulatedSpeciesMagnitude(i));
                Average_response_data=[Average_response_data; ResponseChange];
                
            elseif numel(PeakMagnitude)==0 && numel(TroughMagnitude)==0 % if no peaks or troughs found use maximum value in case response has become saturated
                use=max(RelevantData);
                ResponseChange=use-StimulatedSpeciesMagnitude(i);
                Average_response_data=[Average_response_data; ResponseChange];              
            end

        end
        
        MeanChangeStimulation=mean(Average_response_data);
        PerturbationData=[PerturbationData MeanChangeStimulation]; %extract info for each parameter

        %Return parameter to original value
        modelObj.parameters(p).value=original_parameter;
        
    end
    
    MetDat=[median(PerturbationData) skewness(PerturbationData)];
    MetaData=[MetaData; MetDat];
    PerturbationData=sort(PerturbationData);
    Results=[Results; PerturbationData]; %extract info for each model
    
end


%Plot results
figure(1)
imagesc(Results)
caxis([-1 1])   %colourjet mapping to data
xticks([])
xticklabels([])
yticks([])
yticklabels([])
xlabel('Parameters')
ylabel('Network structures')
set(gca,'FontWeight','bold','fontsize',30)
colorbar
colormap(jet)
     
%Save data and Figures
%filename = 'StimulationChanges.xlsx';
%xlswrite(filename,Results,1);
%h=figure(1);
%savefig(['StimulationChanges'])

