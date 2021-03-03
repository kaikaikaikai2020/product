%Öğ²½À©ÕÅµÄstring
classdef strAdd < handle
    properties
        str1
    end
    methods
        function obj = strAdd(str0)
            if nargin<1
                str0 = [];
            end
            if ~iscell(str0)&&~isempty(str0)
                str0 = {str0};
            end
            [m,n] = size(str0);
            if m<n
                str0 = str0';
            end
            obj.str1 = str0;
        end
        function obj = A(obj,str0)      
            if ~iscell(str0) && ~isempty(str0)
                str0 = {str0};
            end
            obj.str1 = [obj.str1;str0];
        end
    end
	methods(Static)
		function strs = batch_cell_str(strs)
			ind = cellfun(@isnumeric,strs);
            strs(ind) = cellfun(@num2str,strs(ind),'UniformOutput',false);
            strs = strjoin(strs,'\r\n');
		end
	end
end