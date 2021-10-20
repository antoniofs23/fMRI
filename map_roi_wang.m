function [] = map_roi_wang(subjname)

% map wang atlas to locate FEF
%
% ROI for FEF seems better (more accurate when comparing to anatomical landmarks) than Glasser
%
% Wang, L., Mruczek, R. E., Arcaro, M. J., & Kastner, S. (2015). 
% Probabilistic maps of visual topography in human cortex. 
% Cerebral cortex, 25(10), 3911-3931.
%%
% cd to subject
cd(strcat('/Applications/freesurfer/subjects/',subjname));

% copy files to the subject's directory if they dont exist
if ~isfile('lh.Kastner2015.mgz')
    !cp ~/Desktop/Toolboxes/fMRI/lh.Kastner2015.mgz .
    !cp ~/Desktop/Toolboxes/fMRI/rh.Kastner2015.mgz .
end

% if files have already been mapped just open them in freeview
if ~isfile(strcat('FEF_K15_',subjname,'.mgz'))

setenv('MATLAB_SHEL','/bin/bash')

addpath(genpath('/Applications/freesurfer/matlab/'))
setenv('PATH','/usr/local/fsl/bin:/Applications/freesurfer/bin:/usr/bin:/bin')  % this is to tell matlab where freesurfer scripts are
setenv('FREESURFER_HOME','/Applications/freesurfer');  % this is to tell where freesurfer folder is
setenv('SUBJECTS_DIR','/Applications/freesurfer/subjects');  % this is to tell where freesurfer's subject folder is
[~,fsdir] = system('echo $SUBJECTS_DIR'); % find freesurfer directory.
fsdir = fsdir(1:end-1); % there is extra space at the end of echo when called from matlab

freesurfer_init %sometimes needed to run freesurfer scripts in matlab

fr = fileread('ROIfiles_Labeling_Kastner2015.txt'); % read the file that has ROI names with corresponding integers
rois = regexp(fr, '[^\n]*', 'match'); % read text

selected_ROI = 26; % that's FEF ->rois
disp(rois{selected_ROI})

%% Turn into a single file left + right
my_selected_ROI = rois{selected_ROI};

space = strfind(my_selected_ROI,' '); % seperate number from roi name

myROInumber = str2num(my_selected_ROI(1:space-1));
myROIname = my_selected_ROI(space+1:end)

% Let's seperate it from the atlas into a single file for left and right hemisphere

hemi = {'lh';'rh'}

for h = 1 : length(hemi)
    
    %/Applications/freesurfer/matlab/MRIread
    tmp = MRIread(sprintf('%s.Kastner2015.mgz',hemi{h}))
    
    tmp.vol = tmp.vol == 25; % find FEF in the atlas file
    
    MRIwrite(tmp,sprintf('%s.FEF_K15.mgz',hemi{h}))
    
end

%% Let's look at the ROI on the fsaverage space (just run uncommented line that is below), this will block matlab (ctrl + c ) to exit
% !tksurferfv fsaverage lh inflated -overlay rh.FEF_K15.mgz -fminmax 0.1 1

%% Let's map the FEF roi from the fsavearge space to the native space of your subject (I suppose it's Nina)

for h = 1 : length(hemi)
    system(sprintf('%smri_surf2surf --srcsubject fsaverage --trgsubject %s --hemi %s --sval %s.%s.mgz --tval %s.%s_native_K15.mgz' ...
        ,freesurfer_string,subjname,hemi{h},hemi{h},'FEF_K15',hemi{h},'FEF'))
    
end

%% You can look up FEF on the native surface of your subject, this will also block matlab
%system(sprintf('tksurferfv %s lh pial -overlay lh.FEF_native_K15.mgz -fminmax 0.1 1',subjname))

%% I think you need the ROI on the volume so let's map it on the volume space

for h = 1 : length(hemi)
    
    system(sprintf('%smri_surf2vol --hemi %s --template %s/%s/mri/orig.mgz --outvol %s.FEF_native_vol_K15.mgz --surfval %s.FEF_native_K15.mgz --fillribbon --identity %s' ...
        ,freesurfer_string,hemi{h},fsdir,subjname,hemi{h},hemi{h},subjname))
    
end


%% Let's combine hemispheres

lh = MRIread('lh.FEF_native_vol_K15.mgz');
rh = MRIread('rh.FEF_native_vol_K15.mgz');

lh.vol = lh.vol + rh.vol;

MRIwrite(lh,sprintf('%s_%s.mgz','FEF_K15',subjname))

%% Let's convert it to nifti as it's more commonly used format

system(sprintf('mri_convert %s_%s.mgz %s_%s.nii.gz','FEF_K15',subjname,'FEF_K15',subjname))

%% convert T1 to nifti

% path to T1.mgz
T1.mgz = sprintf('/Applications/freesurfer/subjects/%s/mri/T1.mgz',subjname);

% path to T1 nifti file
T1.nii = sprintf('/Applications/freesurfer/subjects/%s/mri/T1.nii.gz',subjname);

% do it
str = sprintf('mri_convert --resample_type nearest -i %s -o %s', T1.mgz, T1.nii);

system(str)

%% Move all files to subject folder

if ~exist(subjname,'dir')
    mkdir(subjname)
    movefile(sprintf('*%s*',myROIname),sprintf('%s/',subjname))
    
else
    
    
    
end

end

%% You can now look it up on the volume (this will block matlab)
%system(sprintf('freeview -v %s/%s/mri/orig.mgz -colormap heat -v FEF_K15_%s.nii.gz  ',fsdir,subjname,subjname))

%% Surface view
right = ' /Applications/freesurfer/subjects/%s/surf/rh.inflated:overlay=/Applications/freesurfer/subjects/%s/rh.FEF_native_K15.mgz';
left  = ' /Applications/freesurfer/subjects/%s/surf/lh.inflated:overlay=/Applications/freesurfer/subjects/%s/lh.FEF_native_K15.mgz';
full  = strcat('freeview -f',right,' ',left,'  -cam azimuth 100 -ss screenshot1.jpg');
system(sprintf(full,subjname,subjname,subjname,subjname));

%% another surface view
%system(sprintf('tksurferfv %s rh pial -overlay rh.FEF_native_K15.mgz -fminmax 0.1 1',subjname))

%% system(sprintf('tksurferfv %s lh inflated -overlay ./bert/lh.FEF_native.mgz -fminmax 0.1 1',subjname)) INFLATED
%% system(sprintf('tksurferfv %s lh pial -overlay ./bert/lh.FEF_native.mgz -fminmax 0.1 1',subjname)) FOLDER
%% system(sprintf('freeview -v %s/%s/mri/orig.mgz -colormap heat -v %s/FEF_%s.nii.gz  ',fsdir,subjname,subjname,subjname))
