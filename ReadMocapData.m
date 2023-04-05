    classdef ReadMocapData
    
    % Version info: this is the latest version of ReadMocapData, as of
    % January 27 2021. Updates will add more features for the different
    % projects that I am involved in. Please, only use this version. 
    % Usage: load data with the class ReadMocapData, execute 
    % data.PreProcessData followed by data.ProcessData for joint angles.
    
    properties(GetAccess = 'protected', SetAccess = 'private') 
        
        Data = [];
        participants = {};
        participanttags = {};
        sides = {};
        
    end 
    
    properties(GetAccess = 'public', SetAccess = 'protected')
        
        Participants = [];
        Warnings = {}; % under construction
        
    end

    properties(GetAccess = 'public', SetAccess = 'public')

        Parameters = [];
        
    end
    
    methods (Access = public)
        
        function obj = ReadMocapData(FullFilePath)
            % read in mocap data, and initiate all data preprocessing and other
            % manipulations
%             close all
            if nargin < 1
                [FileName,DataDirectory] = uigetfile({'*.mat'});
                FullFilePath = fullfile(DataDirectory,FileName);
            end
            try
                obj.Data = mcread(FullFilePath); 
            catch
                obj.Data = load(FullFilePath);
                fieldnames = fields(obj.Data);
                if numel(fieldnames) ~= 1
                    try    
                        obj.Data = obj.Data.mocapdata;
                    catch
                        try
                            obj.Data = obj.Data.data;
                        catch
                            obj.Data = obj.Data.dancer;
                        end
                    end
                else
                    obj.Data = obj.Data.(fieldnames{1});
                end
            end
            try           
                temp0 = load(FullFilePath);
                temp1.S = temp0.mocapdata;
                temp1 = mcreadmat(temp1);
                temp1 = mcfillgaps(temp1, 100);
                obj.Data = temp1;
            catch
                temp = mcfillgaps(obj.Data, 100);
                obj.Data = temp;
            end
%             obj = PreProcessData(obj); 
%             obj = ProcessData(obj); 
            disp([FullFilePath, ' was succesfully loaded.'])
        end

        function obj = PreProcessData(obj)
            % you can set a ROI in the command line, default is all data
            if ~isfield(obj.Parameters,'ROI')
                obj.Parameters.ROI = [1 size(obj.Data.data,1)];
            end
            % cut and resample the data
            if ~isempty(obj.Parameters.ROI)
                data = obj.Data.data(min(obj.Parameters.ROI):max(obj.Parameters.ROI),:);
                obj.Data.data = data;
            end
            % split data for different subjects
            temp = {};
            for ii = 1:numel(obj.Data.markerName)
                if isempty(extractBefore(obj.Data.markerName{ii},'_'))==0
                    temp{end+1} = extractBefore(obj.Data.markerName{ii},'_'); %#ok<*AGROW>
                elseif isempty(extractBefore(obj.Data.markerName{ii},'_'))==1
                    temp{end+1} = 'void';
                end
            end
            
            obj.participants = unique(temp);
            
            for ii = 1:numel(obj.participants)
                if numel(num2str(ii))<2; obj.participanttags{end+1} = ['P0',num2str(ii)];
                elseif numel(num2str(ii))==2; obj.participanttags{end+1} = ['P',num2str(ii)];
                else; disp(['Only up to 9 participants are supported by this script.']);
                end
            end
            
            for ii = 1:numel(obj.participants)
                eval(['obj.Participants.',obj.participanttags{ii},'.Name = obj.participants{ii};']) 
            end
            
            obj.Parameters.ResamplingFactor = 1;
            obj.Parameters.ROI = [];
            if isempty(obj.Parameters.ROI); obj.Parameters.ROI = [1, size(obj.Data.data,1)]; end 
            obj.sides = [{'L'},{'R'}];
            obj.Parameters.Fs = 120;
            obj.Parameters.T = (obj.Parameters.ROI(1):obj.Parameters.ResamplingFactor:obj.Parameters.ROI(2))/obj.Parameters.Fs;

            for ii = 1:numel(obj.participants)
                cc = find(strcmp(temp,obj.participants{ii}));
                markernames = obj.Data.markerName(cc);   %#ok<NASGU>
                cc = ((cc'-1)*3)+1;  
                cc = sort([cc;cc+1;cc+2]); 
                originaldata = obj.Data.data(obj.Parameters.ROI(1):obj.Parameters.ResamplingFactor:obj.Parameters.ROI(2),cc);   %#ok<NASGU>
                eval(['obj.Participants.',obj.participanttags{ii},'.OriginalData.MarkerNames = markernames;'])
                eval(['obj.Participants.',obj.participanttags{ii},'.OriginalData.Data = originaldata;'])
            end
    
            obj = StabilizeData(obj);
            
        end  
        
        function obj = ProcessData(obj)
            % calculate joint angles and body marker data for every participant
            for ii = 1:numel(obj.participants)
%                 eval(['temp = obj.Participants.',obj.participanttags{ii},'.OriginalData.Data;'])
%                 temp = temp*0 + 1;
%                 eval(['obj.Participants.',obj.participanttags{ii},'.ProcessedData.BodyMarkerData = temp;'])
                try
                    eval(['temp = obj.Participants.',obj.participanttags{ii},'.StabilizedData;'])
%                     eval(['temp = obj.',obj.participants{ii},'.StabilizedData;'])
                    for iii = 1:numel(obj.sides)
                        extra = GetMarkerData(obj,temp,obj.sides{iii}); %#ok<NASGU> % just for testing 
                        jointangledata = CalculateJointAngles(obj,temp,obj.sides{iii}); %#ok<NASGU>
                        eval(['obj.Participants.',obj.participanttags{ii},'.ProcessedData.',obj.sides{iii},'.JointAngleData = jointangledata;'])
                        eval(['obj.Participants.',obj.participanttags{ii},'.ProcessedData.',obj.sides{iii},'.Extra = extra;'])
%                         eval(['obj.Participants.',obj.participanttags{ii},'.StabilizedData = [];'])
%                         eval(['obj.Participants.',obj.participanttags{ii},' = rmfield(obj.Participants.',obj.participanttags{ii},',''StabilizedData'');'])
                    end
                catch
                    obj.Warnings{end+1,1} = ['There are not enough data are available to process ',obj.participants{ii}, '!!!'];
                end
            end          
            obj.Data = [];
        end
       
        function obj = PostProcessData(obj)
            for ii = 1:numel(obj.participants)
                try
                    eval(['temp = obj.Participants.',obj.participanttags{ii},'.OriginalData;'])
                    SomeNewParameters = FindParameters(obj,temp); %#ok<NASGU>
                    eval(['obj.Participants.',obj.participanttags{ii},'.PostProcessedData = SomeNewParameters;'])
                catch
                    obj.Warnings{end+1,1} = ['There are not enough data are available to process ',obj.participants{ii}, '!!!'];
                end
            end
            
            for ii = 1:numel(obj.participants)
                try
                    eval(['temp = obj.Participants.',obj.participanttags{ii},'.PostProcessedData;'])
                    subfields = fields(temp);
                    subdiv = obj.AlignSignals(temp.(subfields{1}));
                    eval(['obj.Participants.',obj.participanttags{ii},'.SubDiv = subdiv;']);
                catch
                end
            end
            
            for ii = 1:numel(obj.participants)
                try
                    eval(['temp = obj.Participants.',obj.participanttags{ii},'.PostProcessedData;'])
                    subfields = fields(temp);
                    eval(['subdiv = obj.Participants.',obj.participanttags{ii},'.SubDiv;']);
                    for iii = 1:numel(subfields)
                        eval(['subtemp = obj.Participants.',obj.participanttags{ii},'.PostProcessedData.',subfields{iii},';'])
                        L = length(subtemp(subdiv(1):end));
                        Rows = subdiv(2)-subdiv(1);
                        Columns = floor(L/Rows); 
                        reshapedsignal = reshape(subtemp(subdiv(1):(subdiv(1)+Rows*Columns-1)),Rows,Columns);
                        eval(['obj.Participants.',obj.participanttags{ii},'.PostProcessedData.',subfields{iii},' = [];'])
                        eval(['obj.Participants.',obj.participanttags{ii},'.PostProcessedData.',subfields{iii},'.FID = subtemp;'])
                        eval(['obj.Participants.',obj.participanttags{ii},'.PostProcessedData.',subfields{iii},'.ReshapedSignal = reshapedsignal;'])
                    end
                catch
                end
            end
                    
%                     for ii = 1:numel(obj.participants)
%                 for iii = 1:numel(subfields)
%                     eval(['subtemp = obj.Participants.',obj.participanttags{ii},'.ExtractedFeatures.',subfields{iii},';'])
%                     L = length(subtemp(subdiv(1):end));
%                     Rows = subdiv(2)-subdiv(1);
%                     Columns = floor(L/Rows); 
%                     reshapedsignal = reshape(subtemp(subdiv(1):(subdiv(1)+Rows*Columns-1)),Rows,Columns);
%                     eval(['obj.Participants.',obj.participanttags{ii},'.ExtractedFeatures.',subfields{iii},' = [];'])
%                     eval(['obj.Participants.',obj.participanttags{ii},'.ExtractedFeatures.',subfields{iii},'.FID = subtemp;'])
%                     eval(['obj.Participants.',obj.participanttags{ii},'.ExtractedFeatures.',subfields{iii},'.ReshapedSignal = reshapedsignal;'])
%                 end
%             end 
            
%             for ii = 1:numel(obj.participants)
%                 for iii = 1:numel(subfields)
%                     eval(['temp1 = obj.Participants.',obj.participanttags{1},'.ExtractedFeatures.',subfields{iii},'.ReshapedSignal;'])
%                     localmax = [];
%                     for iiii = 1:size(reshapedsignal,2)
%                         localmax(iiii,1) = obj.FindLags(temp1(:,iiii),temp2(:,iiii));
%                     end
%                     eval(['obj.Participants.SomeResults.',subfields{iii},'.R1 = localmax;'])
%                 end  
%             end
        end
        
        function plotJCS(obj,JointName,tt)
            
            for ii = 1:numel(obj.participants)
                
                eval(['temp = obj.Participants.',obj.participanttags{ii},'.StabilizedData;'])
%                 eval(['temp = obj.',obj.participants{ii},'.StabilizedData;'])
                
                for iii = 1:numel(obj.sides)
                    PlotHelper(obj,obj.participants{ii},temp,tt);
                    side = obj.sides{iii};

                    eval(['[Oout,Xout,Yout,Zout] = JCS_',JointName,'(obj,temp,side);'])

                    [MarkerData] = GetMarkerData(obj,temp,side);

                    Xin = MarkerData.XYZ(:,1:3:end);
                    Yin = MarkerData.XYZ(:,2:3:end);
                    Zin = MarkerData.XYZ(:,3:3:end);

                    plot3([Oout(tt,1), Oout(tt,1)+Xout(tt,1)],[Oout(tt,2), Oout(tt,2)+Xout(tt,2)],[Oout(tt,3), Oout(tt,3)+Xout(tt,3)],'g-o','LineWidth',2)
                    plot3([Oout(tt,1), Oout(tt,1)+Yout(tt,1)],[Oout(tt,2), Oout(tt,2)+Yout(tt,2)],[Oout(tt,3), Oout(tt,3)+Yout(tt,3)],'c-o','LineWidth',2)
                    plot3([Oout(tt,1), Oout(tt,1)+Zout(tt,1)],[Oout(tt,2), Oout(tt,2)+Zout(tt,2)],[Oout(tt,3), Oout(tt,3)+Zout(tt,3)],'b-o','LineWidth',2)
                    plot3(Oout(tt,1),Oout(tt,2),Oout(tt,3),'mo','LineWidth',2)
                    plot3(Xin(tt,:),Yin(tt,:),Zin(tt,:),'.')
                    
                end
            end
        end
        
    end
    
    methods (Access = protected)

        function obj = StabilizeData(obj)
            % stabilize data
            for ii = 1:numel(obj.participants)
                eval(['temp = obj.Participants.',obj.participanttags{ii},'.OriginalData;'])
                
                data_stabilized_R = temp.Data;
                data_stabilized_L = temp.Data;
                [Ot,Xt,Yt,Zt,sentinel] = JCS_thorax(obj,temp,obj.sides{2});    %#ok<ASGLU>
                
                eval(['obj.Participants.',obj.participanttags{ii},'.Sentinel = sentinel;'])
                
                if sum(~isnan(Ot(:))) > 0
                
                    TM = -Ot; % translation, center dataset around origin Ot.

                    for index = 1:size(temp.Data,1)
                        lst1 = [Xt(index,:);...
                            Yt(index,:);...
                            Zt(index,:)];

                        lst2 = [norm(Xt(index,:)) 0 0;...
                            0 norm(Yt(index,:)) 0;...
                            0 0 norm(Zt(index,:))];

                        M=[0 0 0; ...
                            0 0 0; ...
                            0 0 0];

                        outerproduct = zeros(size(lst1,2),size(lst2,2));
                        for index2 = 1:size(lst1,1)
                            x = lst1(index2,:);
                            y = lst2(index2,:);
                            for i1=1:length(x)
                                for i2=1:length(y)
                                    outerproduct(i1,i2) = x(i1)*y(i2);
                                end
                            end
                            M = M + outerproduct;
                        end

                        N11 = M(1,1) + M(2,2) + M(3,3);
                        N22 = M(1,1) - M(2,2) - M(3,3);
                        N33 = -M(1,1) + M(2,2) - M(3,3);
                        N44 = -M(1,1) - M(2,2) + M(3,3);
                        N12 = M(2,3) - M(3,2);
                        N13 = M(3,1) - M(1,3);
                        N14 = M(1,2) - M(2,1);
                        N21 = N12;
                        N23 = M(1,2) + M(2,1);
                        N24 = M(1,3) + M(3,1);
                        N31 = N13;
                        N32 = N23;
                        N34 = M(3,2) + M(2,3);
                        N41 = N14;
                        N42 = N24;
                        N43 = N34;

                        N=[N11,N12,N13,N14; ...
                            N21,N22,N23,N24; ...
                            N31,N32,N33,N34; ...
                            N41,N42,N43,N44];

                        N(isinf(N)|isnan(N)) = 0;

                        [V,~,~] = eig(N);
                        D = eig(N);
                        [~,c]=max(D);

                        temp_ = rotatepoint(quaternion(V(:,c)'),[temp.Data(index,1:3:end)'+TM(index,1),temp.Data(index,2:3:end)'+TM(index,2),temp.Data(index,3:3:end)'+TM(index,3)]);

                        offset = 10000; % cheapo way to make sure everything is in the right quadrant

                        data_stabilized_R(index,1:3:end) = temp_(:,1)' + offset;
                        data_stabilized_R(index,2:3:end) = temp_(:,3)' + offset;
                        data_stabilized_R(index,3:3:end) = temp_(:,2)' + offset; 

                        data_stabilized_L(index,1:3:end) = temp_(:,1)' + offset;
                        data_stabilized_L(index,2:3:end) = -temp_(:,3)' + offset;
                        data_stabilized_L(index,3:3:end) = temp_(:,2)' + offset; 
                    end

                    temp.Data = [];
                    temp.Data(:,:,1) = data_stabilized_R;
                    temp.Data(:,:,2) = data_stabilized_L;
                    eval(['obj.Participants.',obj.participanttags{ii},'.StabilizedData = temp;'])
%                     eval(['obj.',obj.participants{ii},'.StabilizedData = temp;'])

                end
            end
        end 
        
        function [Oout,Xout,Yout,Zout,sentinel] = JCS_thorax(obj,temp,side)

            % ISB definition
            % Ot = incisure jugularis (IJ); not measured
            % Yt = line connecting midpoint processus xiphoideus (PX)-8th thoracic 
            % vertebra (T8) and midpoint IJ-7th cervical vertebra (C7), pointing up.
            % Zt = line perpendicular to plane IJ, C7 and midpoint PX-T8
            % Xt = line perpendicular to plane Yt-Zt
            % IPEM definition, using Qualisys marker setup guide
            % Ot

            MarkerData = GetMarkerData(obj,temp,side);

            Oout = MarkerData.IJ;

            Yt1 = (MarkerData.T8 + MarkerData.PX)/2;
            Yt2 = (Oout + MarkerData.C7)/2;

            % Yt
            Yout = Yt2-Yt1;
            Yout = 100*Yout./vecnorm(Yout,2,2);

            % Zt
            Zout = cross(MarkerData.IJ-Yt1,MarkerData.C7-Yt1); 
            Zout = 100*Zout./vecnorm(Zout,2,2);

            % Xt
            Xout = cross(Yt2-Yt1,Zout);
            Xout = 100*Xout./vecnorm(Xout,2,2);

            sentinel = ones(1,size(Oout,1));
            sentinel(isinf(sum(sum([Oout,Xout,Yout,Zout],2),2))|isnan(sum(sum([Oout,Xout,Yout,Zout],2),2))) = 0;

        end

        function [GH, residuals] = LocateShoulder(obj,temp,side)
            % in order to determine shoulder angles, we need to located the shoulder
            % joint first. I propose that most markers on the upper arm, rotate around
            % one point of rotation: the shoulder joint. If we fix the body according to
            % a reference frame (the JCS of the thorax), we can estimate the location
            % of thE GH on the mocap data. This technique is also reported by Stokdijk
            % et al.

            % close all
            % clc
            MarkerData = GetMarkerData(obj,temp,side);

            HX = MarkerData.H(:,1);
            HY = MarkerData.H(:,2);
            HZ = MarkerData.H(:,3);

            cc = find(~isnan(sum(HX+HY+HZ,2)));

            HX = HX(cc,:);
            HY = HY(cc,:);
            HZ = HZ(cc,:);

            nPoint = size(HX,2);
            nFrame = size(HX,1);
            MPM = [HX(:), HY(:), HZ(:)];
            A = cat(2, 2 .* MPM, zeros(nFrame * nPoint, nPoint));

            for iN = 1:nPoint
              A(((iN - 1) * nFrame + 1):(iN * nFrame), 3 + iN) = 1;
            end

            b = sum(MPM .* MPM, 2);
            x = A \ b;
            GH = transpose(x(1:3));
            residuals = A * x - b;
            residuals = mean(abs(residuals));

        end

        function [Oout,Xout,Yout,Zout] = JCS_humerus(obj,temp,side)

            [MarkerData] = GetMarkerData(obj,temp,side);

            [GH, ~] = LocateShoulder(obj,temp,side); % could also be added as an input
            [~,~,Yf,~] = JCS_forearm(obj,temp,side);

            Oout = repmat(GH,size(MarkerData.US,1),1);
            % Yh2 = (EL + EM)/2 - GH;
            Yout = GH - MarkerData.EL;
            Yout = 100*Yout./vecnorm(Yout,2,2);
            Zout = cross(-Yout,Yf); % sign convention?
            Zout = 100*Zout./vecnorm(Zout,2,2); 
            Xout = cross(Zout,Yout); % sign convention?
            Xout = 100*Xout./vecnorm(Xout,2,2);

        end

        function [Oout,Xout,Yout,Zout] = JCS_forearm(obj,temp,side)

        [MarkerData] = GetMarkerData(obj,temp,side);

            % Yf = (EL + EM)/2 - US;
            Yout = MarkerData.EL - MarkerData.US;
            Yout = 100*Yout./vecnorm(Yout,2,2);
            Oout = MarkerData.US;
            % Xf = cross(US-(EL + EM)/2, US-RS);
            Xout = cross(-MarkerData.US+MarkerData.EL, MarkerData.US-MarkerData.RS); % ***
            Xout = 100*Xout./vecnorm(Xout,2,2);
            Zout = cross(Yout,Xout);
            Zout = 100*Zout./vecnorm(Zout,2,2);

        end

        function [Oout,Xout,Yout,Zout] = JCS_radius(obj,temp,side)

            [MarkerData] = GetMarkerData(obj,temp,side);

            Oout = (MarkerData.US + MarkerData.RS)/2;
            Yout = MarkerData.EL-Oout; % more or less following the ISB
            Yout = 100*Yout./vecnorm(Yout,2,2);
            % there is no sigmoid notch in our marker set! assuming that the 
            % end result is a Zr pointing laterally, in the plane formed by the US, RS
            % and EL. So I decide to define Xr first, as a vector perpendicular to this
            % plane. Zr is perpendicular to Xr and Yr.
            Xout = cross(MarkerData.US-MarkerData.EL,MarkerData.RS-MarkerData.EL);
            Xout = 100*Xout./vecnorm(Xout,2,2);
            Zout = cross(Yout,Xout);
            Zout = 100*Zout./vecnorm(Zout,2,2);

        end

        function [Oout,Xout,Yout,Zout] = JCS_metacarpal(obj,temp,side)

            [MarkerData] = GetMarkerData(obj,temp,side);

            % due to sime missing markers, no fully independent JCS can be determined for
            % metacarpal and radius bones. Thse JCSs will be only useful for
            % flexion/extension and abduction/adduction.

            Oout = MarkerData.MC;
            Yout = MarkerData.RS - Oout; % more or less following the ISB
            Yout = 100*Yout./vecnorm(Yout,2,2);
            % there is no sigmoid notch in our marker set! assuming that the 
            % end result is a Zr pointing laterally, in the plane formed by the US, RS
            % and EL. So I decide to define Xr first, as a vector perpendicular to this
            % plane. Zr is perpendicular to Xr and Yr.
            Xout = cross(Yout, MarkerData.US - Oout);
            Xout = 100*Xout./vecnorm(Xout,2,2);
            Zout = cross(Xout,Yout);
            Zout = 100*Zout./vecnorm(Zout,2,2);

        end

        function [Oout,Xout,Yout,Zout] = JCS_tibialfibula(obj,temp,side)

            [MarkerData] = GetMarkerData(obj,temp,side);
            % Marker choice is not consistent throughout the dataset. Sometimes KneeIn
            % and KneeOut is recorded, sometimes AnkleIn and AnkleOut is recorded. I
            % came up with different JCS based on most common marker choice in the
            % dataset. Also, ISB is not clear in ankle JCS definition. Not ISB-conform!

            % So, according to ISB MM and LM define the torsional plane, and the origin
            % IM. Only LM is present in the data (MM is optional). We replace this with the
            % marker HeelBack. 

            IM = (MarkerData.MM+MarkerData.LM)/2;
            IC = (MarkerData.LC+MarkerData.MCk)/2;
            Oout = IM;
            Yout = IC - IM;
            Zout = cross(Yout, MarkerData.LC - MarkerData.MCk);
            Xout = cross(Yout,-Zout);
            Xout = 100*Xout./vecnorm(Xout,2,2);
            Yout = 100*Yout./vecnorm(Yout,2,2);
            Zout = 100*Zout./vecnorm(Zout,2,2);

        end

        function [Oout,Xout,Yout,Zout] = JCS_calcaneus(obj,temp,side)

            [MarkerData] = GetMarkerData(obj,temp,side);
            % 

            IM = (MarkerData.MM+MarkerData.LM)/2;
            Oout = IM; 
            IT = nanmean(cat(3,MarkerData.FFo,MarkerData.FFi,MarkerData.TT),3);
            Xout = IT - Oout;
            Yout = cross(Xout, MarkerData.LM - MarkerData.MM); 
            Zout = cross(-Xout,Yout);
            Xout = 100*Xout./vecnorm(Xout,2,2);
            Yout = 100*Yout./vecnorm(Yout,2,2);
            Zout = 100*Zout./vecnorm(Zout,2,2);
            
        end

        function [Oout,Xout,Yout,Zout] = JCS_femur(obj,temp,side)

            [MarkerData] = GetMarkerData(obj,temp,side);
            % again problems: not enough markers on the femur to define a full
            % coordinate system. We only define (roughly) the origin: the marker on the
            % knee, and the Y-axis, or the axes along the femur.

            Oout = MarkerData.MCk;
            Yout = MarkerData.T - Oout;
            Yout = 100*Yout./vecnorm(Yout,2,2);
            Xout = nan(size(Yout));
            Zout = nan(size(Yout));

        end

        function [Oout,Xout,Yout,Zout] = JCS_hips(obj,temp,side)

            [MarkerData] = GetMarkerData(obj,temp,side);
            % origin cannot be defined, because this is the head of the femur.

            if strcmpi(side,'R') == 1
                Oout = (MarkerData.ASIS_R + MarkerData.PSIS_R)/2;
                Zout = MarkerData.ASIS_R - MarkerData.ASIS_L;
                MS = MarkerData.ASIS_R - (MarkerData.PSIS_R + MarkerData.PSIS_L)/2;
                proj_MS_Zout = dot(MS,Zout,2)./vecnorm(Zout,2,2) .* Zout./vecnorm(Zout,2,2);
                Xout = MS - proj_MS_Zout;
                Yout = cross(-Zout,Xout);
                Xout = 100*Xout./vecnorm(Xout,2,2);
                Yout = 100*Yout./vecnorm(Yout,2,2);
                Zout = 100*Zout./vecnorm(Zout,2,2);    
            elseif strcmpi(side,'L') == 1
                Oout = (MarkerData.ASIS_L + MarkerData.PSIS_L)/2;
                Zout = MarkerData.ASIS_L - MarkerData.ASIS_R;
                MS = MarkerData.ASIS_L - (MarkerData.PSIS_R + MarkerData.PSIS_L)/2;
                proj_MS_Zout = dot(MS,Zout,2)./vecnorm(Zout,2,2) .* Zout./vecnorm(Zout,2,2);
                Xout = MS - proj_MS_Zout;
                Yout = cross(-Zout,Xout);
                Xout = 100*Xout./vecnorm(Xout,2,2);
                Yout = 100*Yout./vecnorm(Yout,2,2);
                Zout = 100*Zout./vecnorm(Zout,2,2);    
            else
                error('Please indicate if you want to analyze the right side (side argument = ''R'') or the left side (side argument = ''L'')')
            end

        end
        
        function [Oout,Xout,Yout,Zout] = JCS_neck(obj,temp,side)
    
            [MarkerData] = GetMarkerData(obj,temp,side);

            Oout = MarkerData.C7;
            Yout = MarkerData.HT - Oout;
            Zout = MarkerData.HR - MarkerData.HL;
            Xout = cross(Yout,Zout);
            Xout = 100*Xout./vecnorm(Xout,2,2);
            Yout = 100*Yout./vecnorm(Yout,2,2);
            Zout = 100*Zout./vecnorm(Zout,2,2); 

        end
        
        function MarkerData = GetMarkerData(obj,temp,side)
            
            MarkerNames = {'All'}; MarkerData.XYZ = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = ({'Chest'}); MarkerData.IJ = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = ({'WaistLBack', 'WaistRBack'}); T8 = GetMarkerDataColumns(obj,MarkerNames,temp,side,'mean');
            MarkerNames = ({'BackL', 'BackR'}); T8 = cat(3,T8,GetMarkerDataColumns(obj,MarkerNames,temp,side,'mean'));
            MarkerNames = ({'SpineTop','SpineThoracic2'}); T8 = cat(3,T8,GetMarkerDataColumns(obj,MarkerNames,temp,side)); MarkerData.T8 = nanmean(T8,3);
            MarkerNames = ({'Chest','WaistLFront','WaistRFront'}); MarkerData.PX = GetMarkerDataColumns(obj,MarkerNames,temp,side,'nanmean');
            MarkerNames = ({'SpineTop','SpineThoracic2'}); MarkerData.C7 = GetMarkerDataColumns(obj,MarkerNames,temp,side);
%             MarkerNames = ({'LShoulderTop', 'RShoulderTop'}); MarkerData.C7 = GetMarkerDataColumns(obj,MarkerNames,temp,side,'mean'); % this modification is because there is no SpineTop in Shalan's dataset!!! This is not correct!!! warning warning warning
            MarkerNames = strcat(side,{'ElbowOut'}); MarkerData.EL = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'ElbowIn'}); MarkerData.EM = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'WristIn'}); MarkerData.US = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'WristOut'}); MarkerData.RS = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'Arm'}); MarkerData.H = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'HandOut','Hand2'}); MarkerData.MC = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'HeelBack'}); MarkerData.MM = GetMarkerDataColumns(obj,MarkerNames,temp,side); % AnkleIn would be better
            MarkerNames = strcat(side,{'AnkleOut'}); MarkerData.LM = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'Shin','ShinFrontHigh'}); MarkerData.LC = GetMarkerDataColumns(obj,MarkerNames,temp,side); % Should be KneeIn ideally
            MarkerNames = strcat(side,{'KneeOut'}); MarkerData.MCk = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'ForefootIn'}); MarkerData.FFi = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'ForefootOut'}); MarkerData.FFo = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'ToeTip'}); MarkerData.TT = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = strcat(side,{'Thigh','ThighFrontLow'}); MarkerData.T = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = ({'WaistRBack','WaistR'}); MarkerData.PSIS_R = GetMarkerDataColumns(obj,MarkerNames,temp,side); 
            MarkerNames = ({'WaistLBack','WaistL'}); MarkerData.PSIS_L = GetMarkerDataColumns(obj,MarkerNames,temp,side); 
            MarkerNames = ({'WaistRFront'}); MarkerData.ASIS_R = GetMarkerDataColumns(obj,MarkerNames,temp,side); 
            MarkerNames = ({'WaistLFront'}); MarkerData.ASIS_L = GetMarkerDataColumns(obj,MarkerNames,temp,side); 
            MarkerNames = ({'HeadFront'}); MarkerData.HF = GetMarkerDataColumns(obj,MarkerNames,temp,side); 
            MarkerNames = ({'HeadTop'}); MarkerData.HT = GetMarkerDataColumns(obj,MarkerNames,temp,side); 
            MarkerNames = ({'HeadL'}); MarkerData.HL = GetMarkerDataColumns(obj,MarkerNames,temp,side);
            MarkerNames = ({'HeadR'}); MarkerData.HR = GetMarkerDataColumns(obj,MarkerNames,temp,side);

        end
        
        function [XYZ,obj] = GetMarkerDataColumns(obj,MarkerNames,temp,side,flag)

            if strcmpi(side,obj.sides{2}) == 1
                sideArg = 1;
            elseif strcmpi(side,obj.sides{1}) == 1
                sideArg = 2;
            else 
                error('Hm.. something is wrong, buddy..')
            end
            
            if size(temp.Data,3)>1
                X = temp.Data(:,1:3:end,sideArg);
                Y = temp.Data(:,2:3:end,sideArg);
                Z = temp.Data(:,3:3:end,sideArg);
            else
                X = temp.Data(:,1:3:end);
                Y = temp.Data(:,2:3:end);
                Z = temp.Data(:,3:3:end);
            end
            
            cc = [];
            for k = 1:numel(temp.MarkerNames)
                if find(endsWith(temp.MarkerNames(k),MarkerNames)) 
                    cc = cat(1,cc,k);
                end
            end
                    
            if nargin == 4
                if strcmp(MarkerNames,'All') == 0
                    if length(cc)>1
                        obj.Warnings{end+1,1} = ['The naming of body markers is ambigous.']; %#ok<*NBRAK>
                    end
                    if isempty(cc)
                        XYZ = nan(size(temp.Data,1),3); 
                        for index = 1:numel(MarkerNames)
                            obj.Warnings{end+1,1} = ['Marker with label: ',MarkerNames{index}, ' was not found.'];
                        end
                    else
                        XYZ = [X(:,cc),Y(:,cc),Z(:,cc)];
                    end
                elseif strcmp(MarkerNames,'All') == 1
                    if size(temp.Data,3)>1
                        XYZ = temp.Data(:,:,sideArg);
                    else
                        XYZ = temp.Data;
                    end
                end                
            elseif nargin > 4
                if strcmp(flag,'mean') == 1
                    if isempty(cc)
                        XYZ = nan(size(temp.Data,1),3);
                    else
                        XYZ = [mean(X(:,cc),2),mean(Y(:,cc),2),mean(Z(:,cc),2)];
                    end
                elseif strcmp(flag,'nanmean') == 1
                    if isempty(cc)
                        XYZ = nan(size(temp.Data,1),3);
                    else
                        XYZ = [nanmean(X(:,cc),2),nanmean(Y(:,cc),2),nanmean(Z(:,cc),2)];
                    end
                end
            end
        end  
        
        function jointangledata = CalculateJointAngles(obj,temp,side)
            [~,Xt,Yt,Zt,~] = JCS_thorax(obj,temp,side);
            LocateShoulder(obj,temp,side); % find shoulder "joint"
            [~,Xh2,Yh2,Zh2] = JCS_humerus(obj,temp,side); % JCS humerus 
            [~,Xf,Yf,~] = JCS_forearm(obj,temp,side); % JCS forearm
            e1 = Zh2; e1 = e1./vecnorm(e1,2,2); e3 = Yf; e3 = e3./vecnorm(e3,2,2); e2 = cross(Zh2,Yf); e2 = e2./vecnorm(e2,2,2); % define rotational axes
            jointangledata.Elbow.FE = obj.findangle(e2,Yh2,e1,90);  % flexion (+)/hyperextension (-), in order to follow convention
            jointangledata.Elbow.PS = obj.findangle(e2,Xf,-e3,90); % pronation (+)/supination (-), with neutral position knuckles in same plane humerus-elbow-forearm
            [~,~,Yr,Zr] = JCS_radius(obj,temp,side); % JCS radius
            [~,~,Ym,~] = JCS_metacarpal(obj,temp,side); % JCS metacarpal bones (~hand)
            e1 = Zr; e1 = e1./vecnorm(e1,2,2); e3 = Ym; e3 = e3./vecnorm(e3,2,2); e2 = cross(Zr,Ym); e2 = e2./vecnorm(e2,2,2); % define rotational axes
            jointangledata.Wrist.FE = obj.findangle(e2,Yr,-e1,-90); % flexion/extension, flexion towards the palm is positive
            jointangledata.Wrist.AA = obj.findangle(e1,e3,e2,-90); % abduction/adduction, deflection to the ulna (to the outside of the hand with knuckles up) is positive
            e1 = Yt; e1 = e1./vecnorm(e1,2,2); e3 = Yh2; e3 = e3./vecnorm(e3,2,2); e2 = cross(e1,e3); e2 = e2./vecnorm(e2,2,2); % define rotational axes
            jointangledata.Shoulder.AA = obj.findangle(e2,Xt,e1); % abduction (0 deg) <-> forward flexion (90 deg)           
            jointangledata.Shoulder.E = obj.findangle(e1,e3,-e2); % elevation (-)
            jointangledata.Shoulder.PS = obj.findangle(e2,Xh2,-e3); % pronation/supination, to inside (+) outside (-)
            [~,Xtf,Ytf,~] = JCS_tibialfibula(obj,temp,side);
            [~,~,Ycn,Zcn] = JCS_calcaneus(obj,temp,side);
            e1 = Xtf; e1 = e1./vecnorm(e1,2,2); e3 = Ycn; e3 = e3./vecnorm(e3,2,2); e2 = cross(e1,e3); e2 = e2./vecnorm(e2,2,2);
            jointangledata.Ankle.DP = obj.findangle(e1,e3,e2,- 90); % dorsiflexion (+) or plantar- flexion (-)
            jointangledata.Ankle.R = obj.findangle(-e2,Zcn,e3); % internal rotation (+) or external rotation (-)
            jointangledata.Ankle.E = obj.findangle(e2,Ytf,e1,90); % inversion (+) or eversion (-)
            [~,~,Yfm,~] = JCS_femur(obj,temp,side);
            e1 = Yfm; e1 = e1./vecnorm(e1,2,2); e3 = Ytf; e3 = e3./vecnorm(e3,2,2); e2 = cross(e1,e3); e2 = e2./vecnorm(e2,2,2); 
            jointangledata.Knee.FE = obj.findangle(e1,e3,e2); % flexion or extension knee
            [~,~,Yhip,Zhip] = JCS_hips(obj,temp,side);
            e1 = Zhip; e1 = e1./vecnorm(e1,2,2); e3 = Yfm; e3 = e3./vecnorm(e3,2,2); e2 = cross(e1,e3); e2 = e2./vecnorm(e2,2,2); 
            jointangledata.Hip.FE = obj.findangle(-e2,Yhip,e1); % flexion or extension hip
            jointangledata.Hip.AA = obj.findangle(e1,e3,e2,-90); % abduction or adduction hip
            [~,~,Yhd,Zhd] = JCS_neck(obj,temp,side);
            e1 = Zt; e1 = e1./vecnorm(e1,2,2); e3 = Yhd; e3 = e3./vecnorm(e3,2,2); e2 = cross(e1,e3); e2 = e2./vecnorm(e2,2,2); 
            jointangledata.Neck.AA = obj.findangle(e2,Yt,e1,90); % yes
            jointangledata.Neck.BB = obj.findangle(e1,e3,e2,-90); % wobble
            jointangledata.Neck.GG = obj.findangle(e2,Zhd,e3,90); % no
            e1 = Zhip; e1 = e1./vecnorm(e1,2,2); e3 = Yt; e3 = e3./vecnorm(e3,2,2); e2 = cross(e1,e3); e2 = e2./vecnorm(e2,2,2); 
            jointangledata.Spine.AA = obj.findangle(e2,Yhip,e1,90); % bow
            jointangledata.Spine.BB = obj.findangle(e1,e3,e2,-90); % lean left/right
            jointangledata.Spine.GG = obj.findangle(e2,Zt,e3,90); % twist
        end
        
        function PlotHelper(obj,participant,temp,t)
            %PlotHelper Summary of this function goes here
            %   Detailed explanation goes here

%             figure;
            hold on
            title(participant)
%             hold on

            % Head
            MarkerNames = strcat({'HeadTop'});
            HeadTop = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'HeadFront'});
            HeadFront = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'HeadL'});
            HeadL = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'HeadR'});
            HeadR = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            plot3([HeadFront(t,1),HeadL(t,1)],[HeadFront(t,2),HeadL(t,2)],[HeadFront(t,3),HeadL(t,3)],'r-o')
            plot3([HeadL(t,1),HeadR(t,1)],[HeadL(t,2),HeadR(t,2)],[HeadL(t,3),HeadR(t,3)],'r-o')
            plot3([HeadR(t,1),HeadFront(t,1)],[HeadR(t,2),HeadFront(t,2)],[HeadR(t,3),HeadFront(t,3)],'r-o')
            plot3([HeadFront(t,1),HeadTop(t,1)],[HeadFront(t,2),HeadTop(t,2)],[HeadFront(t,3),HeadTop(t,3)],'r-o')
            plot3([HeadL(t,1),HeadTop(t,1)],[HeadL(t,2),HeadTop(t,2)],[HeadL(t,3),HeadTop(t,3)],'r-o')
            plot3([HeadR(t,1),HeadTop(t,1)],[HeadR(t,2),HeadTop(t,2)],[HeadR(t,3),HeadTop(t,3)],'r-o')

            % Waist
            MarkerNames = strcat({'WaistLBack','WaistL'});
            WaistLBack = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'WaistRBack','WaistR'});
            WaistRBack = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'WaistLFront'});
            WaistLFront = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'WaistRFront'});
            WaistRFront = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            plot3([WaistLBack(t,1),WaistRBack(t,1)],[WaistLBack(t,2),WaistRBack(t,2)],[WaistLBack(t,3),WaistRBack(t,3)],'r-o')
            plot3([WaistRBack(t,1),WaistRFront(t,1)],[WaistRBack(t,2),WaistRFront(t,2)],[WaistRBack(t,3),WaistRFront(t,3)],'r-o')
            plot3([WaistRFront(t,1),WaistLFront(t,1)],[WaistRFront(t,2),WaistLFront(t,2)],[WaistRFront(t,3),WaistLFront(t,3)],'r-o')
            plot3([WaistLFront(t,1),WaistLBack(t,1)],[WaistLFront(t,2),WaistLBack(t,2)],[WaistLFront(t,3),WaistLBack(t,3)],'r-o')

            % Back % Chest
            MarkerNames = strcat({'SpineTop','SpineThoracic2'});
            SpineTop = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'Chest'});
            Chest = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'BackL'});
            BackL = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            MarkerNames = strcat({'BackR'});
            BackR = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{2});
            plot3([SpineTop(t,1),BackL(t,1)],[SpineTop(t,2),BackL(t,2)],[SpineTop(t,3),BackL(t,3)],'r-o')
            plot3([BackL(t,1),BackR(t,1)],[BackL(t,2),BackR(t,2)],[BackL(t,3),BackR(t,3)],'r-o')
            plot3([BackR(t,1),SpineTop(t,1)],[BackR(t,2),SpineTop(t,2)],[BackR(t,3),SpineTop(t,3)],'r-o')
            plot3([SpineTop(t,1),Chest(t,1)],[SpineTop(t,2),Chest(t,2)],[SpineTop(t,3),Chest(t,3)],'r-o')
            plot3([BackL(t,1),Chest(t,1)],[BackL(t,2),Chest(t,2)],[BackL(t,3),Chest(t,3)],'r-o')
            plot3([BackR(t,1),Chest(t,1)],[BackR(t,2),Chest(t,2)],[BackR(t,3),Chest(t,3)],'r-o')

            for iii = 1:numel(obj.sides)
                
                % Arm
                MarkerNames = strcat(obj.sides{iii},{'ShoulderTop'});
                ShoulderTop = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'ShoulderBack'});
                ShoulderBack = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});   
                MarkerNames = strcat(obj.sides{iii},{'Arm'});
                Arm = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});   
                MarkerNames = strcat(obj.sides{iii},{'ElbowOut'});
                ElbowOut = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'ElbowIn'});
                ElbowIn = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'WristIn'});
                WristIn = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'WristOut'});
                WristOut = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'HandOut','Hand2'});
                HandOut = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                plot3([ShoulderTop(t,1),ShoulderBack(t,1)],[ShoulderTop(t,2),ShoulderBack(t,2)],[ShoulderTop(t,3),ShoulderBack(t,3)],'r-o')
                plot3([ShoulderBack(t,1),Arm(t,1)],[ShoulderBack(t,2),Arm(t,2)],[ShoulderBack(t,3),Arm(t,3)],'r-o')
                plot3([Arm(t,1),ElbowOut(t,1)],[Arm(t,2),ElbowOut(t,2)],[Arm(t,3),ElbowOut(t,3)],'r-o')
                plot3([ElbowOut(t,1),ElbowIn(t,1)],[ElbowOut(t,2),ElbowIn(t,2)],[ElbowOut(t,3),ElbowIn(t,3)],'r-o')
                plot3([ElbowIn(t,1),ShoulderTop(t,1)],[ElbowIn(t,2),ShoulderTop(t,2)],[ElbowIn(t,3),ShoulderTop(t,3)],'r-o')
                plot3([ElbowOut(t,1),WristOut(t,1)],[ElbowOut(t,2),WristOut(t,2)],[ElbowOut(t,3),WristOut(t,3)],'r-o')
                plot3([ElbowIn(t,1),WristIn(t,1)],[ElbowIn(t,2),WristIn(t,2)],[ElbowIn(t,3),WristIn(t,3)],'r-o')
                plot3([WristOut(t,1),WristIn(t,1)],[WristOut(t,2),WristIn(t,2)],[WristOut(t,3),WristIn(t,3)],'r-o')
                plot3([WristOut(t,1),HandOut(t,1)],[WristOut(t,2),HandOut(t,2)],[WristOut(t,3),HandOut(t,3)],'r-o')

                % Leg
                MarkerNames = strcat(obj.sides{iii},{'Thigh','ThighFrontLow'});
                Thigh = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii}); 
                MarkerNames = strcat(obj.sides{iii},{'KneeOut'});
                KneeOut = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});   
                MarkerNames = strcat(obj.sides{iii},{'KneeIn'});
                KneeIn = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});    
                MarkerNames = strcat(obj.sides{iii},{'Shin','ShinFrontHigh'});
                Shin = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'HeelBack'});
                HeelBack = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'AnkleOut'});    
                AnkleOut = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'ForefootIn'});    
                ForefootIn = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'ForefootOut'});
                ForefootOut = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                MarkerNames = strcat(obj.sides{iii},{'ToeTip'});
                ToeTip = GetMarkerDataColumns(obj,MarkerNames,temp,obj.sides{iii});
                plot3([Thigh(t,1),KneeOut(t,1)],[Thigh(t,2),KneeOut(t,2)],[Thigh(t,3),KneeOut(t,3)],'r-o')
                plot3([KneeOut(t,1),KneeIn(t,1)],[KneeOut(t,2),KneeIn(t,2)],[KneeOut(t,3),KneeIn(t,3)],'r-o')
                plot3([KneeIn(t,1),Thigh(t,1)],[KneeIn(t,2),Thigh(t,2)],[KneeIn(t,3),Thigh(t,3)],'r-o')
                plot3([KneeOut(t,1),Shin(t,1)],[KneeOut(t,2),Shin(t,2)],[KneeOut(t,3),Shin(t,3)],'r-o')
                plot3([KneeIn(t,1),Shin(t,1)],[KneeIn(t,2),Shin(t,2)],[KneeIn(t,3),Shin(t,3)],'r-o')
                plot3([Shin(t,1),AnkleOut(t,1)],[Shin(t,2),AnkleOut(t,2)],[Shin(t,3),AnkleOut(t,3)],'r-o')
                plot3([AnkleOut(t,1),HeelBack(t,1)],[AnkleOut(t,2),HeelBack(t,2)],[AnkleOut(t,3),HeelBack(t,3)],'r-o')
                plot3([AnkleOut(t,1),ForefootOut(t,1)],[AnkleOut(t,2),ForefootOut(t,2)],[AnkleOut(t,3),ForefootOut(t,3)],'r-o')
                plot3([ForefootOut(t,1),ToeTip(t,1)],[ForefootOut(t,2),ToeTip(t,2)],[ForefootOut(t,3),ToeTip(t,3)],'r-o')
                plot3([ToeTip(t,1),ForefootIn(t,1)],[ToeTip(t,2),ForefootIn(t,2)],[ToeTip(t,3),ForefootIn(t,3)],'r-o')

            end

            [XYZ] = GetMarkerDataColumns(obj,'All',temp,obj.sides{2});
            
%             centerpicture = nanmean(cat(3,WaistLBack,WaistRBack,WaistLFront,WaistRFront),3);
            centerpicture = [nanmean(XYZ(:,1:3:end),2),nanmean(XYZ(:,2:3:end),2),nanmean(XYZ(:,3:3:end),2)];
            axis([centerpicture(t,1)-1000, centerpicture(t,1)+1000, ...
                centerpicture(t,2)-1000, centerpicture(t,2)+1000, ...
                centerpicture(t,3)-1000, centerpicture(t,3)+1000])
            view([90 0 0])
            axis square 
            hold off
        
        end    
        
        function SomeNewParameters = FindParameters(obj,temp)
            
            X = temp.Data(:,1:3:end);
            Y = temp.Data(:,2:3:end);
            Z = temp.Data(:,3:3:end);
            
            % deviation center of gravity
            cog = [nanmean(X,2),nanmean(Y,2)];
            for index = 1:size(X,1)
                sel = find(Z(index,:)<40);
                base(index,:) = [nanmean(X(index,sel)); nanmean(Y(index,sel))];
            end
            deviation = sqrt(sum(cog - base,2).^2);
            SomeNewParameters.DeviationGravityCenter = deviation;
            
            % kinetic energy
            v = sqrt(diff(X,[],1).^2 + diff(Y,[],1).^2 + diff(Z,[],1).^2);
            v = nansum(v,2);
            SomeNewParameters.KineticEnergy = v;

            % Z acceleration
            cc = [23,24,25,26]; % waist markers
            a = nansum(diff(Z(:,cc),2,1),2);
            SomeNewParameters.ZAcceleration = a;    
            
            % Z velocity
            v = nansum(diff(Z(:,cc),1,1),2);
            SomeNewParameters.ZVelocity = v;  
            
            % overall stretch
            dcog = sqrt([nansum((X-nanmean(X,2)).^2,2) + ...
                nansum((Y-nanmean(Y,2)).^2,2) + ...
                nansum((Z-nanmean(Z,2)).^2,2)]); 
            SomeNewParameters.Stretch = dcog;  
            
           	% pirouettes
            z = (X(:,cc) - nanmean(X,2)) + 1i*(Y(:,cc) - nanmean(Y,2));
            pangle = real(z(:,1)) + 1i*imag(z(:,1));
            SomeNewParameters.Pirouettes = pangle;  
            
            % bounce
            b = nansum(Z(:,cc),2);
            SomeNewParameters.Bounce = b; 
            
            % z-component
            zc = Z - repmat(nanmean(Z),size(Z,1),1);
            zc = zc./repmat(nanstd(zc),size(zc,1),1);
            SomeNewParameters.ZComponent = zc; 
            
            % retrieve data left and right
            MarkerData_L = GetMarkerData(obj,temp,obj.sides{1});
            MarkerData_R = GetMarkerData(obj,temp,obj.sides{2});
            
            % calculate center of mass
            CenterOfMass = [nanmean(MarkerData_R.XYZ(:,1:3:end),2),nanmean(MarkerData_R.XYZ(:,2:3:end),2),nanmean(MarkerData_R.XYZ(:,3:3:end),2)];
            
            % calculate distance between hands
            SomeNewParameters.Distance_Hands_Mutual_1 = sqrt(nansum((MarkerData_R.MC - MarkerData_L.MC).^2,2));
            SomeNewParameters.Distance_Hands_Mutual_2 = sqrt(nansum((MarkerData_R.US - MarkerData_L.US).^2,2));
            SomeNewParameters.Distance_Hands_Mutual_3 = sqrt(nansum((MarkerData_R.RS - MarkerData_L.RS).^2,2));
            
            % calculate distance between foot tip and center of mass
            SomeNewParameters.Distance_ToeTip_CenterOfMass_R = sqrt((MarkerData_R.FFi(:,1) - CenterOfMass(:,1)).^2 + (MarkerData_R.FFi(:,2) - CenterOfMass(:,2)).^2);
            SomeNewParameters.Distance_ToeTip_CenterOfMass_L = sqrt((MarkerData_R.FFi(:,1) - CenterOfMass(:,1)).^2 + (MarkerData_L.FFi(:,2) - CenterOfMass(:,2)).^2);
            
            % calculate distance between foot tip and floor
            SomeNewParameters.Distance_ToeTip_Floor_R  = abs(MarkerData_R.FFi(:,3));
            SomeNewParameters.Distance_ToeTip_Floor_L  = abs(MarkerData_L.FFi(:,3));
            
        end
    end   
    
    methods (Static, Access = private)
        
        function Output = findangle(A,B,N,offset)
            if nargin == 3
                Output.Angle = rad2deg(atan2(dot(cross(A,B),N,2),dot(A,B,2)));
                Output.AngularVelocity = diff(Output.Angle,1,1);
                Output.AngularAcceleration = diff(Output.Angle,2,1);
            elseif nargin == 4
                Output.Angle = rad2deg(atan2(dot(cross(A,B),N,2),dot(A,B,2))) + offset;
                Output.AngularVelocity = diff(Output.Angle,1,1);
                Output.AngularAcceleration = diff(Output.Angle,2,1);
            end
        end
        
        function Output = AlignSignals(Z)
%             close all
            temp = Z;
            temp(isnan(Z)) = 0;
            XC = (xcorr(temp,'coeff')); llim = 100;
            if length(Z)>3000
                b = XC((length(Z)+llim):(length(Z)+1000+llim));
            else
                b = XC(length(Z):end);
            end
            L = length(Z);
%             plot(b)
            x = 1:length(b);
            A = [x', ones(size(x'))];
            x = A\b;
            XC = b-A*x;
%             XC = highpass(XC,0.005);
            sel1 = find(diff(sign(diff(XC)))<0);
            [~,cc1] = sort(XC(sel1),'descend');
            cc1 = min(cc1(1:3)); firstpeak = sel1(cc1);
            llim = firstpeak + 300; ulim = llim + firstpeak; % 300
            sel2 = sel1(sel1>llim & sel1<ulim);
            [~,cc2] = sort(XC(sel2),'descend');
            cc2 = cc2(1);
            secondpeak = sel2(cc2);
%             subplot(312)
%             hold on
%             plot(XC)
%             plot(firstpeak,XC(firstpeak),'rx')
%             plot(firstpeak*2+llim,XC(firstpeak*2+llim),'kx')
%             plot(secondpeak,XC(secondpeak),'go')
%             axis tight
            DP = secondpeak-firstpeak;
            windowsize = [round(0.75*DP):1:round(1.25*DP)];
            signalshifter = [-round(1.5*DP):1:round(1.5*DP)];
            for ii = 1:length(signalshifter)
                for iii = 1:length(windowsize)
                    numwindows = floor(L/windowsize(iii));
                    leftover = L - windowsize(iii)*numwindows;    
                    shiftedsignal = circshift(Z,signalshifter(ii));
                    test = reshape(shiftedsignal(1:(end-leftover)),windowsize(iii),numwindows);
                    [VM,CM] = max(nanmean(test,2));
                    iterator(ii,iii,1) = VM;
                    iterator(ii,iii,2) = CM;
                end
            end
            field1 = iterator(:,:,1);
            iterator(:,:,1);
            [~,I] = max(field1(:));
            [~,I2] = ind2sub(size(field1),I);
            [~,I3] = min(abs(iterator(:,I2,2)- round(windowsize(I2)/2)));
            subdiv = signalshifter(I3):windowsize(I2):L;
            subdiv = subdiv(subdiv>0);
%             subplot(311)
%             hold on
%             plot(Z)
%             plot(subdiv,Z(subdiv),'ro')
%             xlabel('t (s)')
%             ylabel('X(t)')
%             axis tight
            L = length(Z(subdiv(1):end));
            Rows = subdiv(2)-subdiv(1);
            Columns = floor(L/Rows);
            FID_Z = Z(subdiv(1):(subdiv(1)+Rows*Columns-1));
            FID_Z = FID_Z - min(FID_Z) + 1;
            iterator = [];
            signalshifter = [-round(0.5*Rows):1:round(0.5*Rows)];
            for index = 1:length(signalshifter)
                tempdata = reshape(circshift(FID_Z,signalshifter(index)),Rows,Columns);
                p1 = sum(tempdata(1,:)) + sum(tempdata(Rows,:));
                p2 = sum(sum(tempdata(50:end-50,:)));
                iterator(index) = p2/p1;
            end
            [~,m] = max(iterator);
            Output = subdiv + signalshifter(m);
            Output = Output(Output>0 & Output<length(Z));
%             FID_Z = reshape(Z(Output(1):(Output(1)+Rows*Columns-1)),Rows,Columns);
%             subplot(313)
%             plot(FID_Z)
        end  
        
    end
 
end
