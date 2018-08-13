%        
function terminating_pts = find_skel_ends(input_skeleton_image,varargin)

option = 'not testing';
if ~isempty(varargin)
    for n = 1:1:length(varargin)
        if strcmp(varargin{n},'testing') | strcmp(varargin{n},'not testing')
            option = varargin{n};
        else
            error('Error in input option');
        end
    end
end

% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------
% % % if strcmp(option,'testing')
% % % 	disp('Process the skeleton to find the edges of the island');
% % % end
% ---------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------
%Extract the relative 
[Y,X] = find(input_skeleton_image ~= 0);
skeleton_coord = [X,Y];
terminating_pts = [];
% Relative vectors to the side borders from the current pixel
border_vector = [-1 -1;...
                  0 -1;...
                  1 -1;...
                  1  0;...
                  1  1;...
                  0  1;...
                 -1  1;...
                 -1  0];
             
  disp(strcat('find_skel_ends :',num2str(length(skeleton_coord(:,1)))));
  
for n = 1:1:length(skeleton_coord(:,1))
% %     if strcmp(option,'testing')
% %         disp(strcat('find_skel_ends :',num2str(n),'-of-',num2str(length(skeleton_coord(:,1)))));
% %     end
    %Select the current coordinate to be tested
  
    current_coord = skeleton_coord(n,:);
    %Determine the coordinates of the pixels around this pixel
    border_coord = [border_vector(:,1) + current_coord(1),border_vector(:,2) + current_coord(2)];
    pos = find(ismember(border_coord,skeleton_coord,'rows') == 1);
    if isempty(pos)
        % If this pixel is a island, then it is an edge to a island of 1 pixel
        % Save it.
        terminating_pts = [terminating_pts;current_coord];
    else
        %Default assumption: pixel is an edge unless otherwise stated
        present = 0;
        %Test all the pixels around this current pixel
        for m = 1:1:length(pos);
            %For each surrounding pixel, test if there is a corresponding pixel on the other side
            continunity = [4,5,6] + pos(m);
            g = find(continunity > 8);
            if ~isempty(g)
                continunity(g) = mod(continunity(g),9)+1;
            end
            if any(ismember(continunity',pos,'rows'))
                % If any pixels that are on the opposite side of the current pixel
                % are present then this pixel is not a terminating pixel
                present = 1;
                break;
            end
        end
        % If the current pixel does not fullfill the above condition, then it is an edge
        if present == 0
            terminating_pts = [terminating_pts;current_coord];
        end
    end
end


% ---------------------------------------------------------------------------------------------
% If testing then display the skeleton image with the detected end points
% ---------------------------------------------------------------------------------------------
if strcmp(option,'testing')
   % disp('Display the terminating ends');
% 	figure;
%     %imshow(input_skeleton_image);
%     imagesc(input_skeleton_image);colormap(gray)
% 	hold on;
% 	plot(terminating_pts(:,1),terminating_pts(:,2),'r*');
%     title('Red points indicated detected terminal points');
  
    %xlabel('Test mode for "find skel ends.m" function');
end
