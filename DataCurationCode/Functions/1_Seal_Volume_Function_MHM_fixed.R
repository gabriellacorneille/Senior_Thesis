#Daire Carroll, Gothenburg University, 2023, carrolldaire@gmail.com
#1: Approaching a population level assessment of body size in pinnipeds using drones, an early warning of environmental degradation.
#Automatically estimate length, width, elipsoid volume and complex volume for spatial polygons representing the outlines of harbour seals.

if(!require("sf")){
  install.packages("sf",dependencies = TRUE)
}
if(!require("smoothr")){
  install.packages("smoothr",dependencies = TRUE)
}
if(!require("reshape2")){
  install.packages("reshape2",dependencies = TRUE)
}
if(!require("lwgeom")){
  install.packages("lwgeom",dependencies = TRUE)
}

#########################################

seg_no = 7 #how many segments should the poplygon be broken into?

bu = 4 #the bandwidth value for smoothing - which removes or reduces the impact of seal limbs on measurements

#########################################

apex = function(poly){
  
  #make a matrix of the distance between each coordinate in the polygon
  dist_pts = dist(cbind((st_coordinates(poly)[,1]),(st_coordinates(poly)[,2])),
                  method = "euclidean", diag = TRUE, upper = TRUE
  )
  
  #convert to long format
  dist_pts = melt(as.matrix(dist_pts), varnames = c("row", "col"))
  
  #figure out which points have the longest distance between them
  loc = which(dist_pts$value == max(dist_pts$value)) 
  
  #extract those points
  selection = data.frame(rbind(dist_pts[loc[1],],dist_pts[loc[length(loc)],]))
  
  
  #extract xy coordinates of those two points
  pt1x = st_coordinates(poly)[,1][selection[1,1]]
  pt1y = st_coordinates(poly)[,2][selection[1,1]]
  pt1 = cbind(pt1x,pt1y)
  pt1 = st_point(pt1)
  
  pt2x = st_coordinates(poly)[,1][selection[1,2]]
  pt2y = st_coordinates(poly)[,2][selection[1,2]]
  pt2 = cbind(pt2x,pt2y)
  pt2 = st_point(pt2)
  
  #?repeat part of what you just did? turn them into sf point objects
  pt1 = cbind(pt1x,pt1y)
  pt1 = st_point(pt1)
  pt2 = cbind(pt2x,pt2y)
  pt2 = st_point(pt2)
  
  #return multipoint sf object
  pts = st_multipoint(rbind(pt1,pt2))
  
  return(pts)
  
}	#find the two farthest points in a polygon (poly) returned as a multipoint object (pts)

midpoint = function(sf_lines = NULL){
  g = st_geometry(sf_lines)
  g_mids = lapply(g, function(x){
    coords = as.matrix(x)
    get_mids = function(coords){
      dist = sqrt((diff(coords[, 1])^2 + (diff(coords[, 2]))^2))
      dist_mid = sum(dist)/2
      dist_cum = c(0, cumsum(dist))
      end_index = which(dist_cum > dist_mid)[1]
      start_index = end_index - 1
      start = coords[start_index, ]
      end = coords[end_index, ]
      dist_remaining = dist_mid - dist_cum[start_index]
      mid = start + (end - start) * (dist_remaining/dist[start_index])
      return(mid)
    }
    return(get_mids(coords))
  })
  return(unlist(g_mids))
}	#find the midpoint of a spatial line

#########################################

curved_length_vol = function(p,plt){
  
  pol1 = p 
  
  
  pol2 = smooth(pol1, method = "ksmooth", smoothness = bu) #remove fins and flippers
  #return multipoint object from apex function (two farthest apart points)
  flipper_apex = apex(pol1)
  
  if(plt == TRUE){
    plot(st_geometry(pol1))
    plot(st_geometry(flipper_apex ), add = TRUE, col = "red")
  }
  
  
  #turn the polygon into a linestring
  b_line = st_linestring(st_coordinates(pol1)[,1:2])
  
  #define a vector to sample things by?
  round(nrow(st_coordinates(pol2))/10)
  samp_by = 10
  sample_smooth = seq(round(nrow(st_coordinates(pol2))/samp_by),
                      nrow(st_coordinates(pol2)), by = round(nrow(st_coordinates(pol2))/samp_by))
  
  
  ds = c(NA,NA)
  
  for(i in sample_smooth){
    
    #fill in all the points in the polygon
    h1 = densify(pol1,n = 50)
    
    #take points around the perimenter of the polygon, based on vector established above (num of points/10)
    b_pt  = st_point(st_coordinates(pol2)[i,1:2])
    
    #calc distance between the point and the line?
    b_d = st_distance(b_line,b_pt)
    
    #turn point into sf object and correct coordinate system
    b_pt = st_sfc(b_pt)
    st_crs(b_pt) = st_crs(pol1)
    
    #calc distance between original polygon and the new point?
    b_d2 = st_distance(pol1,b_pt)
    
    ds = rbind(ds,c(b_d,b_d2))
  }
  
  #when the point is off of the original polygon, make the first number negative?
  ds[which(ds[,2]>0),1] = ds[which(ds[,2]>0),1]*-1
  #calculate mean of that
  md = mean(na.omit(ds[,1]))
  
  #use it to draw a buffer around the smoothed polygon to make the smoothed polygon bigger??
  pol2 = st_buffer(pol2,md) #change this to 2
  
  if(plt == TRUE){
    plot(st_geometry(pol2), add = TRUE)
  }
  
  smooth_apex = apex(pol2)
  
  #edited this line from flipper_apex to smooth_apex
  if(plt == TRUE){
    plot(st_geometry(smooth_apex), add = TRUE, col = "red")
  }
  #buff_1 = st_distance(st_point(flipper_apex[,1]),st_point(smooth_apex[,1]), by_element = TRUE)		#get distance between apex points in the smooth and original polygon
  #buff_2 = st_distance(st_point(flipper_apex[,2]),st_point(smooth_apex[,2]), by_element = TRUE)		#get distance between apex points in the smooth and original polygon
  #buff = (buff_1 + buff_2)/2
  buff = st_distance(flipper_apex,smooth_apex)		#get distance between apex points in the smooth and original polygon. 
  # does this return the average difference??
  
  #turn smoothed ploygon into a multiline string
  pol2 = st_cast(pol2, "MULTILINESTRING")
  
  
  pline = st_cast(pol2, "LINESTRING")
  
  cuts = apex(pol2)
  
  if(plt == TRUE){
    plot(st_geometry(cuts), add = TRUE, col = "blue")
  }
  
  
  parts = st_collection_extract(st_split(pline, cuts),"LINESTRING")
  
  ##why is parts splitting it into three?
  # plot(st_geometry(pol2))
  plot(st_geometry(parts[1,]), add = TRUE, col = "green")
  plot(st_geometry(parts[2,]), add = TRUE, col = "yellow")
  plot(st_geometry(parts[3,]), add = TRUE, col = "blue")
  
  
  
  lens = c()
  
  for(i in 1:length(row.names(parts))){
    lens[i] = length(st_coordinates(parts[i,]))	
  }
  
  if(length(lens)>2){
    longsection = which(lens == max(lens))
    shortsection =  which(lens == min(lens))
    midsection = which(lens != max(lens) & lens != min(lens)) 
    
    
    if(length(longsection) == 2){
      
      hs = st_zm(st_multipoint(rbind(st_coordinates(parts[longsection[1],]),st_coordinates(parts[shortsection,]))))
      hl = st_zm(st_multipoint(st_coordinates(parts[longsection[2],])))
      
      if(is.na(match(hs[,1],hl[,1]))[1] == TRUE){
        repp2 = cbind(rev(st_coordinates(parts[longsection[1],])[,1]),rev(st_coordinates(parts[longsection[1],])[,2]))
        repp1 = cbind(rev(st_coordinates(parts[shortsection,])[,1]),rev(st_coordinates(parts[shortsection,])[,2]))
        hs = st_zm(st_multipoint(rbind(repp2,repp1)))
      }
      
      ##I think this is throwing an error because it return TRUE TRUE instead of TRUE
      #I'm not sure if we need to make sure both coordinates don't match, I made it so that if either the x or y doesn't match, you reverse the order
      if(st_coordinates(hs)[1,1] != st_coordinates(hl)[1,1] | 
         st_coordinates(hs)[1,2] != st_coordinates(hl)[1,2]){
        XY = cbind(rev(st_coordinates(hs)[,1]),rev(st_coordinates(hs)[,2]))
        hs = st_multipoint(x = XY, dim = "XY")
      }
      
    }else if(length(shortsection) == 2){
      
      hs = st_zm(st_multipoint(rbind(st_coordinates(parts[shortsection[1],]),st_coordinates(parts[shortsection[2],]))))
      hl = st_zm(st_multipoint(st_coordinates(parts[longsection[1],])))
      
      if(is.na(match(hs[,1],hl[,1]))[1] == TRUE){
        repp2 = cbind(rev(st_coordinates(parts[longsection[1],])[,1]),rev(st_coordinates(parts[longsection[1],])[,2]))
        repp1 = cbind(rev(st_coordinates(parts[shortsection,])[,1]),rev(st_coordinates(parts[shortsection,])[,2]))
        hs = st_zm(st_multipoint(rbind(repp2,repp1)))
      }
      
      
      #same error as above (TRUE TRUE instead of TRUE)
      if(st_coordinates(hs)[1,1] != st_coordinates(hl)[1,1] | 
         st_coordinates(hs)[1,2] != st_coordinates(hl)[1,2]){
        XY = cbind(rev(st_coordinates(hs)[,1]),rev(st_coordinates(hs)[,2]))
        hs = st_multipoint(x = XY, dim = "XY")
      }
      
    }else{
      
      #combine the short and mid section together
      
      #MHM extract first and last point in each section line
      mid_1 <- st_coordinates(parts[midsection,])[1,]
      mid_2 <- st_coordinates(parts[midsection,])[nrow(st_coordinates(parts[midsection,])),]
      
      short_1 <- st_coordinates(parts[shortsection,])[1,]
      short_2 <- st_coordinates(parts[shortsection,])[nrow(st_coordinates(parts[shortsection,])),]
      
      d_1_1 <- dist(rbind(mid_1, short_1), method = "euclidean")
      #should this not happen?? means you should go one, rev the other
      
      d_1_2 <- dist(rbind(mid_1, short_2), method = "euclidean")
      #means you should go short - mid
      
      
      d_2_1 <- dist(rbind(mid_2, short_1), method = "euclidean")
      #means you should go mid-short
      
      d_2_2 <- dist(rbind(mid_2, short_2), method = "euclidean")
      #means you should reverse??
      
      
      if(d_2_1 < d_1_2){
        
        hs = st_zm(st_multipoint(rbind(st_coordinates(parts[midsection,]),st_coordinates(parts[shortsection,]))))
        
        
      }else{
        hs = st_zm(st_multipoint(rbind(st_coordinates(parts[shortsection,]),st_coordinates(parts[midsection,]))))
        

      }
      

      

      hl = st_zm(st_multipoint(st_coordinates(parts[longsection,]))) #perhaps we dont need to worry about the h1 h2 ordering HERE
      
      #if the first point of the short and long arcs don't match
      if(is.na(match(hs[,1],hl[,1]))[1] == TRUE){
        
        #take the coordinates of the mid section and reverse them in order?
        repp2 = cbind(rev(st_coordinates(parts[midsection,])[,1]),rev(st_coordinates(parts[midsection,])[,2]))
        #reverse the order of the short section
        repp1 = cbind(rev(st_coordinates(parts[shortsection,])[,1]),rev(st_coordinates(parts[shortsection,])[,2]))
        
        #combine the reversed short and mid section together
        #idk why this is different than the version above
        hs = st_zm(st_multipoint(rbind(repp2,repp1)))
      }
      if(identical(st_coordinates(hs)[1,1:2] , st_coordinates(hl)[1,1:2]) == FALSE) {
        XY = cbind(rev(st_coordinates(hs)[,1]),rev(st_coordinates(hs)[,2]))
        hs = st_multipoint(x = XY, dim = "XY")
      }
      
    }
    
  }else{
    
    hs = st_zm(st_multipoint(st_coordinates(parts[1,])))
    hl = st_zm(st_multipoint(st_coordinates(parts[2,]))) #perhaps we dont need to worry about the h1 h2 ordering HERE
    
    if(identical(st_coordinates(hs)[1,1:2] , st_coordinates(hl)[1,1:2]) == FALSE){
      XY = cbind(rev(st_coordinates(hs)[,1]),rev(st_coordinates(hs)[,2]))
      hs = st_multipoint(x = XY, dim = "XY")
    }
    
  }  
  
  h1_l = length(st_coordinates(hs)[,1])
  h2_l = length(st_coordinates(hl)[,1])
  
  if(plt == TRUE){
    plot(st_geometry(hs),add = TRUE, col = "red")
    plot(st_geometry(hl),add = TRUE, col = "green")
  }
  
  spine = data.frame(matrix(ncol = 2, nrow = 0))
  
  #take the number of coordinates in the short arc, take a sequence from 1:nlength, by the segment number that was defined at the top of the script
  for(i in 1+round(seq(1,length(st_coordinates(hs)[,1]),by = length(st_coordinates(hs)[,1])/seg_no))){ 
    
    tryCatch({
      #take the ith point in short and long, connect them into a line string
      pts = st_multipoint(rbind(st_coordinates(hs)[i,][1:2],st_coordinates(hl)[i,][1:2]))
      l2 = st_cast(pts, "MULTILINESTRING")
      l2 = st_sfc(l2)
      
      if(plt == TRUE){
        plot(st_geometry(l2), add = TRUE, col = "red")
      }
      
      #find midpoint of each section
      mp = midpoint(l2)
      
      if(plt == TRUE){
        plot(st_geometry(st_point(mp)), add = TRUE, pch = "X")
      }
      
      spine = rbind(spine,  mp)
    }, error=function(e){cat("\n")})
  }
  
  #calc distance between the last X on the spine and the two apexes
  m = rbind(spine[length(spine[,1]),])
  m1 = smooth_apex[1,]
  m2 = smooth_apex[2,]
  d_pts1 = dist(rbind(m,m1), method= "euclidean")
  d_pts2 = dist(rbind(m,m2), method= "euclidean")
  
  #add apexes onto the spine
  if(d_pts2>d_pts1){
    spine = rbind(spine,  smooth_apex[1,]) 
    spine = rbind(smooth_apex[2,],  spine)
  }else{
    spine = rbind(spine,  smooth_apex[2,])
    spine = rbind(smooth_apex[1,],  spine)
  } #HERE
  
  colnames(spine) = c("X","Y")
  
  spine_line = data.matrix(spine)
  spine_line = st_linestring(spine_line)
  len = st_length(spine_line) + 2*buff 
  
  if(plt == TRUE){
    plot(st_geometry(spine_line), add = TRUE, col = "blue")
  }
  
  hs = st_sf(geom=st_geometry(hs))
  hl = st_sf(geom=st_geometry(hl))
  st_crs(hs) = st_crs(pol1)
  st_crs(hl) = st_crs(pol1)
  
  spine_pts = st_as_sf(spine, coords = c("X","Y"))
  st_crs(spine_pts) = st_crs(pol1)
  
  #distance from the perimeter to the spine
  spine_perimiter1 = st_distance(hs,st_as_sf(spine_pts))#, by_element = TRUE)
  spine_perimiter2 = st_distance(hl,st_as_sf(spine_pts))#, by_element = TRUE)
  
  ws = spine_perimiter1 + spine_perimiter2
  
  #width is the max distance between the two sides
  width = max(spine_perimiter1) + max(spine_perimiter2)
  
  
  ###MHM - I don't know what this does
  segs = c(buff,rep(0,seg_no),buff)
  for(i in 1:length(st_geometry(spine_pts))-1){
    p1 = spine_pts[i,]
    p2 = spine_pts[i+1,]
    segs[i] = sum(na.omit(c(segs[i], as.numeric(st_distance(p1,p2)))))
  } 
  segs2 = c(segs[1] + segs[2]/2, rep(0,seg_no-2), segs[seg_no]/2 + segs[seg_no+1]  + segs[seg_no+2])
  
  for(i in 2:seg_no-1){
    segs2[i] = segs[i]/2 + segs[i + 1]/2 
  }
  
  vol_cylinder = 0
  for(i in 1: length(segs2)){
    r = as.numeric(ws[i + 1]/2)
    h = segs2[i]
    vol_cylinder = vol_cylinder + (pi*(r^2)*h)
  }
  
  vol_elipsoid = (4/3)*pi*(width/2)*(width/2)*(len/2)
  
  dims = c(len,width,vol_cylinder,vol_elipsoid,pol1$lat,pol1$lon)
  
  print(paste0("Length (m):",len, "Width (m)", width, sep = " "))
  
  return(dims)
  
} 

