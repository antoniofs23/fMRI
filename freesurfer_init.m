
if isempty(strfind(cd,'media')) == 0
    freesurfer_string = 'source $FREESURFER_HOME/SetUpFreeSurfer.sh;'
else
    freesurfer_string = 'source $FREESURFER_HOME/sources.sh;';
end