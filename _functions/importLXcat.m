% MEJ-ESS: An Enhanced Monte Carlo Simulator for Electron Transport in Gaseous Detectors
%
% This program is a derivative work based on METHES, which is copyrighted by Mohamed Rabie.
%
% Original METHES Copyright and License:
%     Copyright (c) 2015, Mohamed Rabie
%     All rights reserved.
%     
%     Redistribution and use in source and binary forms, with or without
%     modification, are permitted provided that the following conditions are
%     met:
%     - Redistributions of source code must retain the above copyright notice, 
%     this list of conditions and the following disclaimer.
%     - Redistributions in binary form must reproduce the above copyright notice, 
%     this list of conditions and the following disclaimer in the documentation 
%     and/or other materials provided with the distribution
%     - Proper reference is made in publications reporting results obtained 
%     using this software. At present, the preferred way to reference METHES is as follows: 
%     M. Rabie and C.M. Franck, "A Monte Carlo collision code for electron transport in low temperature plasmas", 
%     to be submitted to J. Phys. D: Appl. Phys. (2015) 
%     as well as M. Rabie and C.M. Franck, METHES code, downloaded from www.lxcat.net/download/METHES/, 2015
%         
%     THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%     AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%     IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%     ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%     LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%     CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%     SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%     INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%     CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%     ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%     POSSIBILITY OF SUCH DAMAGE.
%
% MEJ-ESS Modifications Copyright:
%     Copyright (c) 2024-2025, [Asser Mohamed]
%     All rights reserved.
%
%     The modifications and enhancements in MEJ-ESS are licensed under the
%     same terms as the original METHES software.
%
% For more information, please refer to the README.md file in the repository.

classdef importLXcat
    
    properties
        
        % cell array of directories of LXcat files
        dir
        % ionization cross section
        ion
        % attachment cross section
        att
        % elastic cross section
        ela
        % excitation cross section
        exc
        % momentum transfer cross section
        mom
        % effective momentum transfer cross section
        eff
        % threshold energies of ionization
        ionThresh
        % threshold energies of exciations
        excThresh
        % total cross secion
        tot
        % vector giving discret energies
        energy
        % plot (interactive=1) or do not plot data (interactive=0)
        interactive
            
    end
    
    methods
        
        
        % =====================================================================
        %> @brief loads all cross sections from LXcat-files
        %> Parameters: [ion,att,ela,exc,ionThresh,excThresh]
        % =====================================================================
        function obj = importXsections(obj)
            % loads all cross sections from LXcat-files 
            % into arrays ion, att, ela, exc, mom and eff 
            % and threshold energies into arrays ionThresh and excThresh
            
            obj = ionization(obj);
            obj = attachment(obj);
            obj = elastic(obj);
            obj = excitation(obj);
            obj = momentum(obj);
            obj = effective(obj);
            
        end
        
        % =====================================================================
        %> @brief loads ionization cross sections from LXcat-files
        %> Parameters: [ion,ionThresh]
        % =====================================================================
	function obj = ionization(obj)
		for j = 1 : length(obj.dir)
			file = dir([obj.dir{j} '/*.txt']);
			fid = fopen([obj.dir{j} '/' file.name]);

			if fid == -1
				error(['Cannot open file: ' obj.dir{j} '/' file.name]);
			end

			lines = {};
			while ~feof(fid)
				line = strtrim(fgetl(fid));  % Remove leading/trailing whitespace
        line = strrep(line, sprintf('\t'), ' ');  % <-- FIXED TAB REPLACEMENT
				line = strrep(line, '\t', ' ');  % Replace tabs with single space
				lines{end+1} = line;
			end
			fclose(fid);

			% Find all indices with 'IONIZATION'
			ionIdx = [];
			for i = 1:length(lines)
				if strcmp(lines{i}, 'IONIZATION')
					ionIdx(end+1) = i;
				end
			end

			if isempty(ionIdx)
				obj.ion{j} = [];
			else
				for k = 1:length(ionIdx)
					idx = ionIdx(k);

					% Read ionization threshold energy (2 lines after 'IONIZATION')
					ionThreshLine = str2double(lines{idx + 2});
					obj.ionThresh{j}{k} = ionThreshLine;

					% Find the next '---' separator after that
					dashStart = -1;
					for i = idx+3:length(lines)
					 if length(lines{i}) >= 3 && strcmp(lines{i}(1:3), '---')
							dashStart = i;
							break;
					 end
					end

					% Find the next '---' after the data block
					dashEnd = -1;
					for i = dashStart+1:length(lines)
						if length(lines{i}) >= 3 && strcmp(lines{i}(1:3), '---')
							dashEnd = i;
							break;
						end
					end

					% Extract cross section data
					sigma = [];
					for i = dashStart+1 : dashEnd-1
						parts = sscanf(lines{i}, '%f');
						if length(parts) >= 2
							sigma(end+1,1:2) = parts(1:2)';
						end
					end
          disp(sigma)
					obj.ion{j}{k} = sigma;
				end
			end
		end
	end

        
        % =====================================================================
        %> @brief loads attachment cross sections from LXcat-files
        %> Parameters: [att]
        % =====================================================================
		function obj = attachment(obj)
		    % loads multiple attachment cross sections from LXcat-files into array att 
		
		    for j = 1 : length(obj.dir)
		
		        file = dir([obj.dir{j} '/*.txt']);
		        fid = fopen([obj.dir{j} '/' file.name]);
		
		        lines = {};
		        tline = fgetl(fid);
		        while ischar(tline)
		            lines{end+1} = tline;
		            tline = fgetl(fid);
		        end
		        fclose(fid);
		
		        att_indices = [];
		        for i = 1:length(lines)
		            if length(lines{i}) >= 10 && strcmp(lines{i}(1:10), 'ATTACHMENT')
		                att_indices(end+1) = i;
		            end
		        end
		
		        if isempty(att_indices)
		            obj.att{j} = [];
		        else
		            k = 1;
		            for a = 1:length(att_indices)
		                att_start = att_indices(a);
		
		                % Find two dashed lines after this ATTACHMENT section
		                first_dash = 0;
		                second_dash = 0;
		                dash_count = 0;
		                for i = att_start:length(lines)
		                    if length(lines{i}) >= 2 && lines{i}(1) == '-' && lines{i}(2) == '-'
		                        dash_count = dash_count + 1;
		                        if dash_count == 1
		                            first_dash = i;
		                        elseif dash_count == 2
		                            second_dash = i;
		                            break;
		                        end
		                    end
		                end
		
		                if first_dash > 0 && second_dash > first_dash
		                    data = [];
		                    for i = first_dash+1 : second_dash-1
		                        parts = sscanf(lines{i}, '%f %f');
		                        if length(parts) == 2
		                            data(end+1, :) = parts';
		                        end
		                    end
		                    obj.att{j}{k} = data;
		                    k = k + 1;
		                end
		            end
		        end
		    end
		end

        
        % =====================================================================
        %> @brief loads elastic cross sections from LXcat-files
        %> Parameters: [ela]
        % =====================================================================
             function obj = elastic(obj)
                % loads elastic cross sections from LXcat-files
                % into array ela

                for j = 1 : length(obj.dir)

                        file = dir([obj.dir{j} '/*.txt']);
                        fid = fopen([obj.dir{j} '/' file.name]);

                        C = textscan(fid, '%s', 'Delimiter', '\n');
                        lines = C{1};

                        found = false;
                        data_started = false;
                        sigma = [];

                        for i = 1:length(lines)
                                line = strtrim(lines{i});

                                if ~found
                                        if length(line) >= 7
                                                if strcmp(line(1:7), 'ELASTIC')
                                                        found = true;
                                                end
                                        end
                                elseif ~data_started
                                        if length(line) >= 3
                                                dash_count = 0;
                                                for c = 1:length(line)
                                                        if line(c) == '-'
                                                                dash_count = dash_count + 1;
                                                        end
                                                end
                                                if dash_count >= 3
                                                        data_started = true;
                                                        continue;
                                                end
                                        end
                                else
                                        % check if we reached end (empty or another dashed line)
                                        if isempty(line)
                                                break;
                                        end

                                        dash_count = 0;
                                        for c = 1:length(line)
                                                if line(c) == '-'
                                                        dash_count = dash_count + 1;
                                                end
                                        end
                                        if dash_count >= 3
                                                break;
                                        end

                                        % read the numeric data
                                        tokens = sscanf(line, '%f %f');
                                        if length(tokens) == 2
                                                sigma(end+1, 1) = tokens(1);
                                                sigma(end, 2) = tokens(2);
                                                
                                        end
                                end
                        end

                        if isempty(sigma)
                                obj.ela{j} = [];
                        else
                                obj.ela{j}{1} = sigma;
                        end

                        fclose(fid);
                end
        end


        

        
        % =====================================================================
        %> @brief loads excitation cross sections from LXcat-files
        %> Parameters: [exc,excThresh]
        % =====================================================================   
		function obj = excitation(obj)
		    % loads excitation cross sections from LXcat-files 
		    % into array exc and excThresh threshold energy into excThresh
		    
		    for j = 1 : length(obj.dir)
		    
		        file = dir([obj.dir{j} '/*.txt']);
		        fid = fopen([obj.dir{j} '/' file.name]);
		    
		        lines = {};
		        tline = fgetl(fid);
		        while ischar(tline)
		            lines{end+1} = tline;
		            tline = fgetl(fid);
		        end
		        fclose(fid);
		    
		        exc_indices = [];
		        for i = 1:length(lines)
		            if length(lines{i}) >= 10 && strcmp(lines{i}(1:10), 'EXCITATION')
		                exc_indices(end+1) = i;
		            end
		        end
		    
		        if isempty(exc_indices)
		            obj.exc{j} = [];
		            obj.excThresh{j} = [];
		        else
		            k = 1;
		            for a = 1:length(exc_indices)
		                exc_start = exc_indices(a);
		                thresh_index = exc_start + 2;
		                thresh_energy = str2num(lines{thresh_index});
		    
		                % Find two dashed lines after this EXCITATION section
		                first_dash = 0;
		                second_dash = 0;
		                dash_count = 0;
		                for i = exc_start:length(lines)
		                    if length(lines{i}) >= 2 && lines{i}(1) == '-' && lines{i}(2) == '-'
		                        dash_count = dash_count + 1;
		                        if dash_count == 1
		                            first_dash = i;
		                        elseif dash_count == 2
		                            second_dash = i;
		                            break;
		                        end
		                    end
		                end
		    
		                if first_dash > 0 && second_dash > first_dash
		                    data = [];
		                    for i = first_dash+1 : second_dash-1
		                        parts = sscanf(lines{i}, '%f %f');
		                        if length(parts) == 2
		                            data(end+1, :) = parts';
		                        end
		                    end
		                    obj.exc{j}{k} = data;
		                    obj.excThresh{j}{k} = thresh_energy;
		                    k = k + 1;
		                end
		            end
		        end
		    end
		end

        % =====================================================================
        %> @brief loads momentum transfer cross sections from LXcat-files
        %> Parameters: [ela]
        % =====================================================================
        function obj = momentum(obj)
            % loads momentum transfer  cross sections from LXcat-files 
            % into array mom
            for j = 1 : length(obj.dir)
                
                file = dir([obj.dir{j} '/*.txt']);
                fid = fopen([obj.dir{j} '/' file.name]);
                
                C = textscan(fid, '%s %s %s %s %s %s %s %s' );
                
                x = strfind(C{1,1}, 'MOMENTUM');
                K = sum(cell2mat(x));
                
                if K == 0
                    
                    obj.mom{j} = [];;
                    
                else
                    
                    
                    next = 1 ;
                    for k = 1 : K
                        sigma = 0;
                        i = next ;
                        while  isempty(x{i})
                            i = i+1;
                        end
                        ind = i;
                        
                        x = strfind(C{1,1}, '---');
                        i = ind;
                        while isempty(x{i})
                            i = i+1;
                        end
                        first = i+1;
                        
                        i = first;
                        while isempty(x{i})
                            i = i+1;
                        end
                        last = i-1;
                        
                        N = length(first:last);
                        for i = 1 : N
                            ind = first+i-1;
                            sigma(i,1) = str2num(cell2mat(C{1,1}(ind)));
                            sigma(i,2) = str2num(cell2mat(C{1,2}(ind)));
                            
                        end
                        
                        obj.mom{j}{k} = sigma;
                        next = last + 2 ;
                    end
                    
                end
                
            end
            
            %%
            
            fclose(fid);
            
        end
        
        % =====================================================================
        %> @brief loads effective momentum transfer cross sections from LXcat-files
        %> Parameters: [ela]
        % =====================================================================
        function obj = effective(obj)
            % loads effective cross sections from LXcat-files 
            % into array eff
            
            for j = 1 : length(obj.dir)
                
                file = dir([obj.dir{j} '/*.txt']);
                fid = fopen([obj.dir{j} '/' file.name]);
                
                C = textscan(fid, '%s %s %s %s %s %s %s %s' );
                
                x = strfind(C{1,1}, 'EFFECTIVE');
                K = sum(cell2mat(x));
                
                if K == 0
                    
                   obj.eff{j} = [];
                    
                else
                    
                    
                    next = 1 ;
                    for k = 1 : K
                        sigma = 0;
                        i = next ;
                        while  isempty(x{i})
                            i = i+1;
                        end
                        ind = i;
                        
                        x = strfind(C{1,1}, '---');
                        i = ind;
                        while isempty(x{i})
                            i = i+1;
                        end
                        first = i+1;
                        
                        i = first;
                        while isempty(x{i})
                            i = i+1;
                        end
                        last = i-1;
                        
                        N = length(first:last);
                        for i = 1 : N
                            ind = first+i-1;
                            sigma(i,1) = str2num(cell2mat(C{1,1}(ind)));
                            sigma(i,2) = str2num(cell2mat(C{1,2}(ind)));
                            
                        end
                        
                        obj.eff{j}{k} = sigma;
                        next = last + 2 ;
                    end
                    
                end
                
            end
            
            %%
            
            fclose(fid);
            
        end
        
        % =====================================================================
        %> @brief converts effective mometum transfer cross section to  
        % elastic cross section
        %> Parameters: [ela]
        % ===================================================================== 
        function obj = effective2elastic(obj)
            % converts effective mometum transfer cross section eff to  
            % elastic cross section ela
            
            % gas species
            for j = 1 : length(obj.dir)
                
                 if isempty(obj.ela{j}) & ~isempty(obj.eff{j})
                     
                % collision types
                for k = 1 : length(obj.eff{j})
                   
                        
                        energy = obj.eff{j}{k}(:,1);
                        sigma = obj.eff{j}{k}(:,2);
                        
                        % formula from Vahedi et al (1995) p.190
                        beta = 1/2*( energy.*log(1+energy) )./...
                            ( energy - log(1+energy) );
                        % for energy = 0 --> beta = 1:
                         beta(energy<1e-6) = 1;
                        
                        %sigma = beta.*sigma;
                        
                        % substract excitation cross sections
                        for l = 1 : length(obj.exc{j})
                            
                            energy_eff = obj.exc{j}{l}(:,1);
                            sigma_eff = obj.exc{j}{l}(:,2);
                            energy_eff = [energy_eff ; 1e6];
                            sigma_eff = [sigma_eff ; sigma_eff(end)];
                            sigma = sigma - interp1(energy_eff,sigma_eff,energy,'linear');
                            
                        end
                        
                        % substract ionization cross sections
                        for l = 1 : length(obj.ion{j})
                            
                            energy_eff = obj.ion{j}{l}(:,1);
                            sigma_eff = obj.ion{j}{l}(:,2);
                            energy_eff = [energy_eff ; 1e6];
                            sigma_eff = [sigma_eff ; sigma_eff(end)];
                            sigma = sigma - interp1(energy_eff,sigma_eff,energy,'linear');
                            
                        end
                        
                        obj.ela{j}{k}(:,1) = energy;
                        obj.ela{j}{k}(:,2) = max(0,sigma);
                        
                    end
                    
                end
                
            end
            
        end
             
         % =====================================================================
        %> @brief Converts momentum transfer cross section to  
        % elastic cross section
        %> Parameters: [ela]
        % ===================================================================== 
        function obj = momentum2elastic(obj)
            % converts momentum transfer mom cross section to  
            % elastic cross section ela
            
            % gas species
            for j = 1 : length(obj.dir)
                
                if isempty(obj.ela{j}) & ~isempty(obj.mom{j})
                    
                % collision types
                for k = 1 : length(obj.mom{j})
                    
                                           
                        energy = obj.mom{j}{k}(:,1);
                        sigma = obj.mom{j}{k}(:,2);
                        
                        % formula from Vahedi et al (1995) p.190
                        beta = 1/2*( energy.*log(1+energy) )./...
                            ( energy -log(1+energy) );
                        % for energy = 0 --> beta = 1:
                        beta(energy<1e-6) = 1;
                        
                        sigma = beta.*sigma;
                        
                        obj.ela{j}{k}(:,1) = energy;
                        obj.ela{j}{k}(:,2) = sigma;
                        
                    end
                    
                end
                
            end
            
        end
           
        % =====================================================================
        %> @brief creates over all energy from LXcat data
        %> Parameters: [energy]
        % ===================================================================== 
        function obj = getEnergy(obj)
            % creates over all vector energy from LXcat data
            
            % gas species
            for j = 1 : length(obj.dir)
                
                % collision types
                sigma = obj.att;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                
                sigma = obj.ela;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                
                sigma = obj.exc;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                
                sigma = obj.ion;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                
                sigma = obj.mom;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                    
                
                % find unique values and sort them
                obj.energy = unique(obj.energy);
                obj.energy = sort(obj.energy);
                
                
            end
            
        end
        
        % =====================================================================
        %> @brief Prepare data for fit by adding data points 
        %> at energy = 0 and max(energy) 
        %> Parameters: [ion,att,exc]
        % ===================================================================== 
        function obj = prepareForFit1(obj)
            % prepare data for fit by adding data points at energy = 0 
            % and max(energy) for arrays ion, att and exc
            
            E_max = max(obj.energy);
             
            obj.ion  = obj.addDataPoints(obj.ion ,E_max);
            obj.att  = obj.addDataPoints(obj.att ,E_max);
            obj.exc  = obj.addDataPoints(obj.exc ,E_max);  

        end
        
        % =====================================================================
        %> @brief prepare data for fit by adding data points 
        %> at energy = 0 and max(energy) 
        %> Parameters: [ela,mom,eff]
        % ===================================================================== 
        function obj = prepareForFit2(obj)
            % prepare data for fit by adding data points at energy = 0 
            % and max(energy) for arrays ela, mom and eff
           
            E_max = max(obj.energy);
            
            if ~isempty(obj.ela)
            obj.ela  = obj.addDataPoints(obj.ela,E_max);
            end
            
            if ~isempty(obj.mom)
            obj.mom  = obj.addDataPoints(obj.mom,E_max );
            end
            
            if ~isempty(obj.eff)
            obj.eff  = obj.addDataPoints(obj.eff,E_max );
            end
            
        end
        
        % =====================================================================
        %> @brief adds data point to cross section sigma at energy = 0 and
        % zero-data point at energy, if sigma is empty 
        % ===================================================================== 
        function [sigma] = addDataPoints(obj,sigma,energy)
            % adds data point to cross section sigma at energy = 0 
                  
            % gas species
            for j = 1 : length(obj.dir)

                if isempty(sigma{j})
                    sigma{j} = {[energy'  0*energy']};
                end
                
                % collision types
                for k = 1 : length(sigma{j});
                    
                    xfit = sigma{j}{k}(:,1);
                    yfit = sigma{j}{k}(:,2);
                    
                    % add addiitonal point at beginning:
                    if xfit(1) > 0
                        xfit = [0 ; xfit];
                        yfit = [yfit(1) ; yfit];
                    end
                     
                    sigma{j}{k} = [xfit  yfit] ;
                    
                end
                
            end
            
        end
        
        % =====================================================================
        %> @brief  sets ionization and excitation threshold energies to zero if 
        % not available  
        %> Parameters: [ionThresh,excThresh]
        % ===================================================================== 
        function obj = fillThresholds(obj)
            % sets ionization and excitation threshold energies ionThresh 
            % and excThresh to zero if not available  
            
            % gas species
            for j = 1 : length(obj.dir)

                % ionization energies:
                if isempty(obj.ion{j})
                    obj.ionThresh{j} = {0};
                end
                
                % ionization energies:
                if isempty(obj.exc{j})
                    obj.excThresh{j} = {0};
                end
                                
            end
            
        end
        
        % =====================================================================
        %> @brief  makes total cross section and corresponding energy 
        % for each species
        %> Parameters: [energy,tot]
        % ===================================================================== 
        function obj = totalXsection(obj)
            % makes total cross section tot and corresponding energy 
            % for each species
           
             for j = 1 : length(obj.dir)
                 
             sigma_ion = sumXsection(obj,obj.ion{j},obj.energy);
             sigma_att = sumXsection(obj,obj.att{j},obj.energy);
             sigma_ela = sumXsection(obj,obj.ela{j},obj.energy);
             sigma_exc = sumXsection(obj,obj.exc{j},obj.energy);
             
             obj.tot{j} = sigma_ion + sigma_att + sigma_ela + sigma_exc;
             
             end
 
        end
        
        % =====================================================================
        %> @brief makes sum of one type of cross section over one species
        % ===================================================================== 
       function sigma_sum = sumXsection(obj, sigma, energy)
    % Makes sum of one type of cross section over one species

    if isempty(sigma)
        sigma_sum = 0;
    else
        sigma_sum = 0;
        for k = 1:length(sigma)
            x = sigma{k}(:,1);
            y = sigma{k}(:,2);

            % Ensure x is unique by removing duplicates
            [x, idx] = unique(x); % Keeps the original order
            y = y(idx); % Update y to match the unique x values

            % If the energy range exceeds x, extend x and y
            if energy(end) > x(end)
                x = [x; energy(end)];
                y = [y; y(end)];
            end

            % Perform interpolation
            sigma_add = interp1(x, y, energy);
            if ~isnan(sigma_add)
                sigma_sum = sigma_sum + sigma_add;
            end
        end
    end
end
        % =====================================================================
        %> @brief plots fitted cross sections from LXcat
        % ===================================================================== 
       function plotXsections(obj, textsize, xscale, yscale)
                % Plots fitted cross sections from LXcat

                if obj.interactive == 1
                        % gas species
                        for j = 1:length(obj.dir)

                                [fig, ax] = obj.GetFigure(xscale, yscale);
                                hold on;
                                grid on;
                                Legend = {};

                                for i = 1:length(obj.ion{j})
                                        plot(obj.ion{j}{i}(:,1), obj.ion{j}{i}(:,2), 'r.-', 'LineWidth', 2);
                                        Legend{end+1} = 'ionization';
                                end

                                for i = 1:length(obj.att{j})
                                        plot(obj.att{j}{i}(:,1), obj.att{j}{i}(:,2), 'b.-', 'LineWidth', 2);
                                        Legend{end+1} = 'attachment';
                                end

                                for i = 1:length(obj.exc{j})
                                        plot(obj.exc{j}{i}(:,1), obj.exc{j}{i}(:,2), 'k.-', 'LineWidth', 2);
                                        obj.text2function(obj.exc{j}{i}(:,1), obj.exc{j}{i}(:,2), [num2str(obj.excThresh{j}{i}) ' eV'], 12);
                                        Legend{end+1} = ['excitation ' num2str(obj.excThresh{j}{i}) ' eV'];
                                end

                                for i = 1:length(obj.ela{j})
                                        plot(obj.ela{j}{i}(:,1), obj.ela{j}{i}(:,2), 'm.', 'LineWidth', 2);
                                        Legend{end+1} = 'elastic';
                                end

                                for i = 1:length(obj.mom{j})
                                        if sum(obj.mom{j}{i}(:,2)) > 0
                                                plot(obj.mom{j}{i}(:,1), obj.mom{j}{i}(:,2), 'g-', 'LineWidth', 2);
                                                Legend{end+1} = 'momentum transfer';
                                        end
                                end

                                for i = 1:length(obj.eff{j})
                                        if sum(obj.eff{j}{i}(:,2)) > 0
                                                plot(obj.eff{j}{i}(:,1), obj.eff{j}{i}(:,2), 'c.-', 'LineWidth', 2);
                                                Legend{end+1} = 'effective momentum transfer';
                                        end
                                end

                                plot(obj.energy, obj.tot{j}, '--', 'LineWidth', 2);
                                Legend{end+1} = 'total';

                                xlabel('energy [eV]', 'fontsize', textsize);
                                ylabel('\sigma  [eV]', 'fontsize', textsize);
                                set(gca, 'fontsize', textsize);
                                legend(Legend, 'location', 'northeast');
                                title(['gas ' num2str(j)]);

                                disp(['Cross sections gas ' num2str(j) ' ok']);
                                disp('Press any key to continue...');
                                pause;
                        end
                end
        end

        % =========================================================================
        %> @brief get figure properties
        % =========================================================================
        function [figure1, axes1] = GetFigure(obj, xscale, yscale)
                % Get figure properties

                % Set default scales: linear
                if nargin < 3
                        xscale = 'linear';
                        yscale = 'linear';
                end

                fontsize = 14;
                linewidth = 0.5;

                figure1 = figure('PaperSize', [29.68 20.98]);
                axes1 = axes('Parent', figure1, ...
                             'XScale', xscale, 'XMinorTick', 'off', ...
                             'YScale', yscale, 'YMinorTick', 'on', ...
                             'LineWidth', linewidth, 'FontSize', fontsize);

                box on;
                grid on;
                xlabel('time [ns]', 'LineWidth', linewidth, 'FontSize', fontsize);
                ylabel('amplitude [A]', 'LineWidth', linewidth, 'FontSize', fontsize);
        end

        % =====================================================================
        %> @brief writes string to x- and y-data in Figure 
        % ===================================================================== 
        function text2function(obj,x,y,string,textsize)
            % writes text "string" on the maximum of function y(x)
            
            y_max = max(y);
            y_max = y_max(1);
            x_max = x(y==max(y));
            x_max = x_max(1);
            text(x_max,y_max, string ,'fontsize',textsize)
            
        end
        
    end
    
end


