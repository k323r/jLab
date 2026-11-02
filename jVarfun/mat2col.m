function[varargout]=mat2col(varargin)
%MAT2COL  Compress NAN-padded matrix data into long columns.
%
%   C=MAT2COL(M), where the columns of matrix M are data vectors of
%   nonuniform length with NANs filling out the gaps in each column,
%   appends together all the columns into a vector C in which isolated
%   NANs mark the transitions between columns.
%
%   [C1,C2,...,CN]=MAT2COL(M1,M2,...,MN) also works for multiple input 
%   matrices of the same size.  In this case the locations of NANs in 
%   M1 are used as the key for appending all the MN into columns.  
%  
%   MAT2COL, COL2MAT, and COLBREAKS together form a system for moving
%   data with segments of nonuniform length rapidly back and forth
%   between a column format and a padded-matrix format.
%
%   MAT2COL(M1,M2,...,MN); with no output arguments overwrites the
%   original input variables.
%
%   See also COL2MAT, COLBREAKS, COL2CELL, CELL2COL.
%   _________________________________________________________________
%   This is part of JLAB --- type 'help jlab' for more information
%   (C) 2000--2016 J.M. Lilly --- type 'help jlab_license' for details        
  
%no loops!

if strcmpi(varargin{1},'--t')
    mat2col_test,return
end

if isempty(varargin{1})
   for i=1:nargout
       varargout{i}=[];
   end
   return
end

i=0;
while i<nargin
	i=i+1;
	mat=varargin{i};

	if i==1
		bpad=~all(isnan(mat(end,:)));
	end
	if bpad
		mat(size(mat,1)+1,:)=nan*ones(size(mat(1,:)));
	end
	mat=mat(:);
	if i==1

		%new matrix should include all data points
		bool=false(size(mat));
		index=find(~isnan(mat));
		bool(index)=ones(size(index));
		index2=find(diff([bool;0])==-1)+1;
		if ~isempty(index2)
			bool(index2)=ones(size(index2));
		end
		index=find(bool);
	end
	mat=mat(index);
	col=mat;
	varargout{i}=col;
end


if nargout>nargin
        varargout{end+1}=index;
	%eval(['c' int2str(nargout) '=index;']);
end

if nargout==0
   %for i=1:nargin
   %   eval(['varargout{',int2str(i),'}=c',int2str(i),';'])
   %end  
  eval(to_overwrite(nargin));
end

function[]=mat2col_test

[c1,c2]=mat2col([1 2;3 nan;5 nan],[10 20;30 nan;50 nan]);
reporttest('MAT2COL multiple inputs, last row not all NaN',aresame([c1 c2],[1 10;3 30;5 50;nan nan;2 20;nan nan]))
