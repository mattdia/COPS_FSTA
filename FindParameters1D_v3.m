function [Parameters] = FindParameters1D_v3(parameters_path)

% upgraded for single loop (SmartScan)
% upgraded for Zurich Instrument Lock-In data  - Gael 5/6/2014

AxisNumber =0;
AxisInitialPosition =0;
AxisStepSize= 0;
AxisNumberOfSteps = 0;
AxisWait = 0;
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

Parameters = zeros(3,1);

parFilename= parameters_path;
fid=fopen(parFilename);
        
     for(i=1:1:20)
	 
         tline = fgetl(fid);
         
         	if( findstr(tline,'Axis number') ~= 0)
			 pos = tline;
			 index = strfind(pos,'r');
			 pos_str_val = pos(1,index+1:end);
			 AxisNumber = str2num(pos_str_val); 
        end
	 
		if( findstr(tline,'Stage initial position [microns]') ~= 0)
			 pos = tline;
			 index = strfind(pos,']');
			 pos_str_val = pos(1,index+1:end);
			 AxisInitialPosition = str2num(pos_str_val); 
        end
        
		
         if( findstr(tline,'step size [microns]') ~= 0)
 			 step = tline;
 			 index = strfind(step,']');
 			 step_str_val = step(1,index+1:end);
 			 AxisStepSize = str2num(step_str_val); 
         end
         
          
        
           if( findstr(tline,'Number of steps') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ps');
  			 Num_str_val = Num(1,index+2:end);
  			 AxisNumberOfSteps = str2num(Num_str_val);
           end
           
         
      
           
            if( findstr(tline,'Wait time [ms]') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,']');
  			 Num_str_val = Num(1,index+1:end);
  			 AxisWait = str2num(Num_str_val); 
            end
           
                       
           
            if( findstr(tline,'average pnts') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'ts');
  			 Num_str_val = Num(1,index+2:end);
  			 avgpts = str2num(Num_str_val); 
            end
           
            if( findstr(tline,'Demod 1 TC') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'TC');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod1 = str2num(Num_str_val); 
            end
            
             if( findstr(tline,'Demod 2 TC') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'TC');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod2 = str2num(Num_str_val); 
             end
            
             if( findstr(tline,'Demod 3 TC') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'TC');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod3 = str2num(Num_str_val); 
             end
            
             if( findstr(tline,'Demod 4 TC') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'TC');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod4 = str2num(Num_str_val); 
             end
            
               if( findstr(tline,'Demod 5 TC') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'TC');
  			 Num_str_val = Num(1,index+2:end);
  			 TC_Demod5 = str2num(Num_str_val); 
               end
         
              if( findstr(tline,'Demod 6 TC') ~= 0)
  			 Num = tline;
  			 index = strfind(Num,'TC');
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

           
     end  
     Parameters(1,1) = AxisNumber;
     
     Parameters(2,1) =  AxisInitialPosition;
   
     Parameters(3,1) = AxisStepSize;
 
     Parameters(4,1) =  AxisNumberOfSteps;
  
     Parameters(5,1) = AxisWait;
  
     Parameters(6,1) =  avgpts;
     Parameters(7,1) = TC_Demod1;
     Parameters(8,1) = TC_Demod2;
     Parameters(9,1) = TC_Demod3;
     Parameters(10,1) = TC_Demod4;
     Parameters(11,1) = TC_Demod5;
     Parameters(12,1) = TC_Demod6;
     Parameters(13,1) = Order_Demod1;
     Parameters(14,1) = Order_Demod2;
     Parameters(15,1) = Order_Demod3;
     Parameters(16,1) = Order_Demod4;
     Parameters(17,1) = Order_Demod5;
     Parameters(18,1) = Order_Demod6;
     
%      Parameters(2,1) = StepSize;
%      Parameters(3,1) = NumberOfSteps;
     fclose(fid);
end