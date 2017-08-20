% Load an image containing the checkerboard pattern

% imageFileName = fullfile(toolboxdir('vision'),'visiondata', 'calibration', 'webcam', 'image4.tif')
  imageFileName =fullfile('C:\','Users','devan','Downloads','varsity','DIC selfstudy','assessmentImages','calimg01.tif')
  I = imread(imageFileName);
 
  % Detect the checkerboard points
  [imagePoints, boardSize] = detectCheckerboardPoints(I);
 
  % Display detected points
  J = insertText(I, imagePoints, 1:size(imagePoints, 1));
  J = insertMarker(J, imagePoints, 'o', 'Color', 'red', 'Size', 5);
  imshow(J);
  title(sprintf('Detected a %d x %d Checkerboard', boardSize));