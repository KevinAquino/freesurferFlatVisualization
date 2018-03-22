function [surfaceHandle fh flatPatch parcel] = visualizeFlatMapFreesurfer(subject,annotationFile,flatPatch,hemi,labelCol)
if(nargin<5)
	txtCol = 'white';
	txtSize = 16;
	txtWeight = 'bold';
else
	txtCol = labelCol.col;
	txtSize = labelCol.siz;
	txtWeight = labelCol.weight;
end

% Load up all the surfaces and flat patches
subjectsFolder = ['/Applications/freesurfer/subjects/'];

flatPatch = read_patch([subjectsFolder subject '/surf/' flatPatch]);
[vertFull facesFull] = read_surf([subjectsFolder subject '/surf/' hemi '.white']);
curvData = read_curv([subjectsFolder subject '/surf/' hemi '.curv']);
[verts label ctab] = read_annotation(annotationFile);

[h h2] = ismember(facesFull,flatPatch.ind);
sum3 = sum(h,2);
fac2 = h2(sum3==3,:);

flatVerts = [flatPatch.x;flatPatch.y].';

fh = figure;
surfaceHandle = patch('Vertices',flatVerts,'Faces',fac2,'CData',curvData(flatPatch.ind+1),'FaceColor','interp','EdgeColor','none','FaceAlpha',0.5);
% surfaceHandle = patch('Vertices',flatVerts,'Faces',fac2,'CData',curvData(flatPatch.ind+1),'FaceColor','interp','EdgeColor','none','FaceAlpha',1);
axis image;
set(fh,'Color',[0 0 0]);


msh_curvature           = -curvData.';
mod_depth               = 0.5;
curvatureColorValues    = ((2*msh_curvature>0) - 1) * mod_depth * 128 + 127.5;
curvData = [curvatureColorValues;curvatureColorValues;curvatureColorValues].';
curvData = curvData(flatPatch.ind+1,:)/255;
set(surfaceHandle,'FaceVertexCdata',curvData);

% % going back and forth between the full and small flattened mesh.
fullToSub = sparse(flatPatch.ind+1,1:flatPatch.npts,1,length(vertFull),length(flatPatch.ind)); 
subToFull = sparse(1:flatPatch.npts,flatPatch.ind+1,1,length(flatPatch.ind),length(vertFull));


% Go to every ROI and get the points that form the ROI and get the area! (skip the first one as it is the unknown region)
for nr=2:181,
	roi = verts(label==ctab.table(nr,5));
	% After we have the ROI, then grab the sub-face list (i.e. all the faces that belong to the mesh inside the ROI)
	[h h2] = ismember(facesFull,roi);
	sum3 = sum(h,2);
	sfac=h2(sum3==3,:);
	% Now get all the edges for the triangles inside the ROI, and sort them
	edg1 = sort(sfac(:,1:2),2);
	edg2 = sort(sfac(:,2:3),2);
	edg3 = sort(sfac(:,[1 3]),2);
	% we want to find all the instances where an edge only belongs to one face, these are the perimeter edges. To do this use setxor i.e. find where edges are not in both
	perim = setxor(setxor(edg1,edg2,'rows'),edg3,'rows');
	% Now to get the proper points, we find all the vertices in there
	pts = roi(unique(perim(:)));
% 	keyboard
	% Now get these vertices described in the full mesh and find where they lie in the submesh by making use of the transformation I described above..
	[nold submeshVertices] = find(fullToSub(pts+1,:));
	% Also get the actual ROI in terms of the flat patch -- this is for use of getting the center of mass for labelling
	[nold flatInds] = find(fullToSub(roi+1,:));
	% now store the actual points as well as the color that is used throughout the maps..
	parcel{nr-1,1} = submeshVertices;	
	parcel{nr-1,2} = ctab.table(nr,1:3);	
	parcel{nr-1,3} = ctab.struct_names{nr}(3:end-4);

	% Here is a fancy trick to find out if we need to have two different labels because of the splitting that may occur when cutting the surface for flattening
	coos = [flatPatch.x(parcel{nr-1,1});flatPatch.y(parcel{nr-1,1})];
    meanPoint=mean(coos,2);
    centeredPoints = coos - (repmat(meanPoint,[1 size(coos,2)]));
   	[~,s,v] = svd(centeredPoints);
	newpts = s*(v.');
	xnewpts = newpts(1,:);   
	rangepts = max(xnewpts)-min(xnewpts);

	% Here is an arbitrary threshold to decide if the ROI is too wide, it is currently set at 80mm, but you may want to change this
	if(rangepts>80)			
        
		% Now do the SVD on the entire ROI if the perimeter's x-dimension in the principal direction is too long!
        [~,roiOnFlat] = find(fullToSub(roi+1,:));
        coos = [flatPatch.x(roiOnFlat);flatPatch.y(roiOnFlat)];
        meanPoint=mean(coos,2);
        centeredPoints = coos - (repmat(meanPoint,[1 size(coos,2)]));

        [~,s,v] = svd(centeredPoints);
        newpts = s*(v.');
        xnewpts = newpts(1,:);
        
		part1 = find(xnewpts<= rangepts/4 + min(xnewpts));
		part2 = find(xnewpts> rangepts*3/4 + min(xnewpts));	

        parcel{nr-1,4} = [mean([flatPatch.x(roiOnFlat(part1));flatPatch.y(roiOnFlat(part1))],2) mean([flatPatch.x(roiOnFlat(part2));flatPatch.y(roiOnFlat(part2))],2)];
	else
		parcel{nr-1,4} = mean([flatPatch.x(flatInds);flatPatch.y(flatInds)],2);
	end	
end


for nr=1:180;
	text_h = 0;
	hold on;
	if(~isempty(parcel{nr}))
		plot(flatPatch.x(parcel{nr,1}),flatPatch.y(parcel{nr,1}),'.','Color',parcel{nr,2}/255,'MarkerSize',6);			
		if(size(parcel{nr,4},2) == 2)						
			text_h(1) = text(parcel{nr,4}(1,1),parcel{nr,4}(2,1),parcel{nr,3},'fontSize',txtSize,'interpreter','none','Color',txtCol,'HorizontalAlignment','center','FontWeight',txtWeight);
			text_h(2) = text(parcel{nr,4}(1,2),parcel{nr,4}(2,2),parcel{nr,3},'fontSize',txtSize,'interpreter','none','Color',txtCol,'HorizontalAlignment','center','FontWeight',txtWeight);
		else
			text_h(1) = text(parcel{nr,4}(1),parcel{nr,4}(2),parcel{nr,3},'fontSize',txtSize,'interpreter','none','Color',txtCol,'HorizontalAlignment','center','FontWeight',txtWeight);
			text_h(2) = text_h(1);
		end
	end
	parcel{nr,5} = text_h;
end

axis off;

