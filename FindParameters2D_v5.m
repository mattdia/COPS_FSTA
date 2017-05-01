function [Parameters] = FindParameters2D_v5(parameters_path)

% upgraded for single loop (SmartScan)
% upgraded for Zurich Instrument Lock-In data  - Gael 5/6/2014
%Upgraded for 6D scans. 8/7/2014

Axis1InitialPosition =0;
Axis2InitialPosition =0;
Axis3InitialPosition =0;
Axis1StepSize= 0;
Axis2StepSize= 0;
Axis3StepSize= 0;
Axis1NumberOfSteps = 0;
Axis2NumberOfSteps = 0;
Axis3NumberOfSteps = 0;
Axis1Wait = 0;
CutOffIndex=0;
TotalNumberOfSteps = 0;
avgpts = 0;
TC_Demod1 = 0;
TC_Demod2 = 0;
TC_Demod3 = 0;
TC_Demod4 = 0;
TC_Demod5 = 0;
TC_Demod6 = 0;
Order_Demod1 = 0;
Order_Demod2 = 0;
Order_Demod3 = 0;
Order_Demod4 = 0;
Order_Demod5 = 0;
Order_Demod6 = 0;
V_init = 0;
V_StepSize=0;
VNumberOfSteps = 0;
aux_init = 0;
aux_StepSize = 0;
aux_NumberOfSteps = 0;
aux2_init = 0;
aux2_StepSize = 0;
aux2_NumberOfSteps = 0; %%%%
Parameters = zeros(34,1);

parFilename= parameters_path;
fid=fopen(parFilename);
        
     for(i=1:1:35)
	 
         tline = fgetl(fid);
	 
		if( findstr(tline,'Axis 1 initial position [microns]') ~= 0)
			 pos = tline;
			 index = strfind(pos,']');
			 pos_str_val = pos(1,index+1:end);
			 Axis1InitialPosition = str2num(pos_str_val); %rate in microns
        end
        
		if( findstr(tline,'Axis 2 initial position [microns]') ~= 0)
			 pos = tline;
			 index = strfind(pos,']');
			 pos_str_val = pos(1,index+1:end);
			 Axis2InitialPosition = str2num(pos_str_val); %rate in microns
        end
        
		if( findstr(tline,'Axis 3 initial position [microns]') ~= 0)
			 pos = tline;
			 index = strfind(pos,']');
			 pos_str_val = pos(1,index+1:end);
			 Axis3InitialPosition = str2num(pos_str_val); %rate in microns
        end
        
         if( findstr(tline,'Axis 1 step size [microns]') ~= 0)
 			 step = tline;
 			 index = strfind(step,']');
 			 step_str_val = step(1,index+1:end);
 			 Axis1StepSize = str2num(step_str_val); %rate in microns
         end
         
          if( findstr(tline,'Axis 2 step size [microns]') ~= 0)
 			 step = tline;
 			 index = strfind(step,']');
 			 step_str_val = step(1,index+1:end);
 			 Axis2StepSize = str2num(step_str_val); %rate in microns
          end
          
          if( findstr(tline,'Axis 3 step size [microns]') ~= 0)
 			 step = tline;
 			 index = strfind(step,']');
 			 step_str_val = step(1,index+1:end);
 			 Axis3StepSize = str2num(step_str_val); %rate in microns
         end
         
        
           if( findstr(tline,'Axis 1 Number of steps') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ps');
  			 Num_str_val = Num(1,index+2:end);
  			 Axis1NumberOfSteps = str2num(Num_str_val); %rate in microns
           end
           
           if( findstr(tline,'Axis 2 Number of steps') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ps');
  			 Num_str_val = Num(1,index+2:end);
  			 Axis2NumberOfSteps = str2num(Num_str_val); %rate in microns
           end
           
           if( findstr(tline,'Axis 3 Number of steps') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ps');
  			 Num_str_val = Num(1,index+2:end);
  			 Axis3NumberOfSteps = str2num(Num_str_val); %rate in microns
           end
           
            if( findstr(tline,'Axis 1 Wait time [ms]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,']');
  			 Num_str_val = Num(1,index+1:end);
  			 Axis1Wait = str2num(Num_str_val); 
            end
           
             if( findstr(tline,'Smart Scan cutoff index') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ex');
  			 Num_str_val = Num(1,index+2:end);
  			 CutOffIndex = str2num(Num_str_val); 
            end
             
           if( findstr(tline,'Number of scan points') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ts');
  			 Num_str_val = Num(1,index+2:end);
  			 TotalNumberOfSteps = str2num(Num_str_val); 
           end
           
            if( findstr(tline,'average pnts') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ts');
  			 Num_str_val = Num(1,index+2:end);
  			 avgpts = str2num(Num_str_val); 
            end
           
            if( findstr(tline,'Demod 1 TC [s]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod1 = str2num(Num_str_val); 
            end
            
             if( findstr(tline,'Demod 2 TC [s]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod2 = str2num(Num_str_val); 
             end
            
             if( findstr(tline,'Demod 3 TC [s]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod3 = str2num(Num_str_val); 
             end
            
             if( findstr(tline,'Demod 4 TC [s]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod4 = str2num(Num_str_val); 
             end
            
               if( findstr(tline,'Demod 5 TC [s]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod5 = str2num(Num_str_val); 
               end
         
              if( findstr(tline,'Demod 6 TC [s]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod6 = str2num(Num_str_val); 
              end
            
              
               if( findstr(tline,'Demod 1 filter slope') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'pe');
  			 Num_str_val = Num(1,index+2:end);
  			 Order_Demod1 = str2num(Num_str_val); 
               end
              
                
               if( findstr(tline,'Demod 2 filter slope') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'pe');
  			 Num_str_val = Num(1,index+2:end);
  			 Order_Demod2 = str2num(Num_str_val); 
               end
              
                
               if( findstr(tline,'Demod 3 filter slope') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'pe');
  			 Num_str_val = Num(1,index+2:end);
  			 Order_Demod3 = str2num(Num_str_val); 
               end
              
                
               if( findstr(tline,'Demod 4 filter slope') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'pe');
  			 Num_str_val = Num(1,index+2:end);
  			 Order_Demod4 = str2num(Num_str_val); 
               end
              
                
               if( findstr(tline,'Demod 5 filter slope') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'pe');
  			 Num_str_val = Num(1,index+2:end);
  			 Order_Demod5 = str2num(Num_str_val); 
               end
              
                
               if( findstr(tline,'Demod 6 filter slope') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'pe');
  			 Num_str_val = Num(1,index+2:end);
  			 Order_Demod6 = str2num(Num_str_val); 
              end

              if( findstr(tline,'V initial [V]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'V]');
  			 Num_str_val = Num(1,index+2:end);
  			 V_init = str2num(Num_str_val); 
               end
 
               if( findstr(tline,'V step size  [V]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'V]');
  			 Num_str_val = Num(1,index+2:end);
  			 V_StepSize = str2num(Num_str_val); 
               end
              
               if( findstr(tline,'V number of steps') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ps');
  			 Num_str_val = Num(1,index+2:end);
  			 VNumberOfSteps = str2num(Num_str_val); 
               end
              
         if( findstr(tline,'Aux initial [microns]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 aux_init = str2num(Num_str_val); 
         end
         
         if( findstr(tline,'Aux step size [microns]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 aux_StepSize = str2num(Num_str_val); 
         end
         
         if( findstr(tline,'Aux number of steps') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ps');
  			 Num_str_val = Num(1,index+2:end);
  			 aux_NumberOfSteps = str2num(Num_str_val); 
         end
         
         if( findstr(tline,'Aux2 initial [microns]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 aux2_init = str2num(Num_str_val); 
         end
         
         if( findstr(tline,'Aux2 step size [microns]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'s]');
  			 Num_str_val = Num(1,index+2:end);
  			 aux2_StepSize = str2num(Num_str_val); 
         end
         
         if( findstr(tline,'Aux2 number of steps') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ps');
  			 Num_str_val = Num(1,index+2:end);
  			 aux2_NumberOfSteps = str2num(Num_str_val); 
         end

     end  
     Parameters(1,1) =  Axis1InitialPosition;
     Parameters(2,1) =  Axis2InitialPosition;
     Parameters(3,1) =  Axis3InitialPosition;
     Parameters(4,1) = Axis1StepSize;
     Parameters(5,1) = Axis2StepSize;
     Parameters(6,1) =  Axis3StepSize;
     Parameters(7,1) =  Axis1NumberOfSteps;
     Parameters(8,1) =  Axis2NumberOfSteps;
     Parameters(9,1) = Axis3NumberOfSteps;
     Parameters(10,1) = Axis1Wait;
     Parameters(11,1) = CutOffIndex;
     Parameters(12,1) = TotalNumberOfSteps;
     Parameters(13,1) =  avgpts;
     Parameters(14,1) = TC_Demod1;
     Parameters(15,1) = TC_Demod2;
     Parameters(16,1) = TC_Demod3;
     Parameters(17,1) = TC_Demod4;
     Parameters(18,1) = TC_Demod5;
     Parameters(19,1) = TC_Demod6;
     Parameters(20,1) = Order_Demod1;
     Parameters(21,1) = Order_Demod2;
     Parameters(22,1) = Order_Demod3;
     Parameters(23,1) = Order_Demod4;
     Parameters(24,1) = Order_Demod5;
     Parameters(25,1) = Order_Demod6;
     Parameters(26,1) = V_init;
     Parameters(27,1) = V_StepSize;
     Parameters(28,1) = VNumberOfSteps;
     Parameters(29,1) = aux_init;
     Parameters(30,1) = aux_StepSize;
     Parameters(31,1) = aux_NumberOfSteps;
     Parameters(32,1) = aux2_init;
     Parameters(33,1) = aux2_StepSize;
     Parameters(34,1) = aux2_NumberOfSteps;
%      Parameters(2,1) = StepSize;
%      Parameters(3,1) = NumberOfSteps;
     fclose(fid);
end