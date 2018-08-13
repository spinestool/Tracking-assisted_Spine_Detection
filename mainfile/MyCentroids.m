 function centroids = MyCentroids(numberOfBlobs,blobMeasurements)
 
 for k = 1:numberOfBlobs
    blobCentroid_vector(k,:)= blobMeasurements(k).Centroid;
 end
 centroids = blobCentroid_vector;
 end
 
