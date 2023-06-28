/************************************************************************
 * 
 *                            Remove cloud 
 *  
 * *********************************************************************/

function bitwiseExtract(value, fromBit, toBit) {
  if (toBit === undefined) toBit = fromBit
  var maskSize = ee.Number(1).add(toBit).subtract(fromBit)
  var mask = ee.Number(1).leftShift(maskSize).subtract(1)
  return value.rightShift(fromBit).bitwiseAnd(mask)
}

function cloudfree_mod09a1(image){
  var qa = image.select('StateQA')
  var cloudState = bitwiseExtract(qa, 0, 1) 
  var cloudShadowState = bitwiseExtract(qa, 2)
  var cirrusState = bitwiseExtract(qa, 8, 9)
  var mask = cloudState.eq(0) // Clear
  .and(cloudShadowState.eq(0)) // No cloud shadow
  .and(cirrusState.eq(0)) // No cirrus
  return image.updateMask(mask)  
}

/************************************************************************
 * 
 *                            data
 * 
 * *********************************************************************/
var roi=table
var year='2021'
var st_date=year+'-06-01'
var en_date=year+'-10-01'


var scale_my=500


//======================================> filter date
// var one_year=true  //  true false
// if (one_year){
//   var filter_data=ee.Filter.date(st_date,en_date)
// }else{
//   var filter_data=ee.Filter.date('2018-01-01','2020-01-01')
// }
var filter=ee.Filter.and(
                ee.Filter.date(st_date,en_date),
                ee.Filter.bounds(roi)
             )

//======================================> imageCollection
var col_a1 = ee.ImageCollection('MODIS/006/MOD09A1').filter(filter)  // wet ndvi ndbsi
var col_lst = ee.ImageCollection('MODIS/006/MOD11A1').filter(filter) // lst
print('col_a1',col_a1)
print('col_lst',col_lst)

//======================================> vision RGB
var trueColor =col_a1.select(['sur_refl_b01', 'sur_refl_b04', 'sur_refl_b03']).median().clip(roi);
var trueColorVis = {
  min: -100.0,
  max: 3000.0,
};
Map.centerObject(roi, 6);
Map.addLayer(trueColor, trueColorVis, 'True Color');



/************************************************************************
 * 
 *                            index
 * 
 * Ref. https://www.mdpi.com/2072-4292/11/20/2345
 * *********************************************************************/
  /*https://www.agriculturejournals.cz/publicFiles/121_2018-JFS.pdf
 sur_refl_b01 red   //620-670nm
 sur_refl_b02 nir1  //841-876nm
 sur_refl_b03 blue  //459-479nm
 sur_refl_b04 green //545-565nm
 sur_refl_b05 nir2  //1230-1250nm
 sur_refl_b06 swir1 //1628-1652nm
 sur_refl_b07 swir2 //2105-2155nm
 */
col_lst=col_lst.map(function(img){
  var new_img=img.select('LST_Day_1km').multiply(0.02).subtract(273.5).clip(roi)
  return new_img.copyProperties(img,img.propertyNames())
})

// calculate index
var a1_index={
   // mNDWI
   mNDWI:function(img){
      var mndwi = img.normalizedDifference(["green","swir1"]);  
      return img.addBands(mndwi.rename("mNDWI"));
   },
   // NDVI
   NDVI:function(img){
      var ndvi = img.normalizedDifference(["nir1","red"]);  //0.85 - 0.88 µm
      return img.addBands(ndvi.rename("NDVI"));
   },
   // WET
   WET:function WET(img){
        var wet =img.expression(
          "float(0.1147*b1+ 0.2489*b2+ 0.2408*b3+ 0.3132*b4- 0.3122*b5- 0.6416*b6- 0.5087*b7)",
        {
          "b1":img.select("red"),
          "b2":img.select("nir1"),
          "b3":img.select("blue"),
          "b4":img.select("green"),
          "b5":img.select("nir2"),
          "b6":img.select("swir1"),
          "b7":img.select("swir2"),
        }
        );
        return img.addBands(wet.multiply(0.0001).rename("WET"));
    },
    // NDBSI
    NDBSI:function NDBSI(img){ // NDBSI
      var ibi = img.expression(
        '((2*swir/(swir+nir))-((nir/(nir+red))+(green/(green+swir))))/((2*swir/(swir+nir))+((nir/(nir+red))+(green/(green+swir))))',
        {
          'swir':img.select('swir1'),//1.57 - 1.65 µm
          'nir':img.select('nir1'),
          'red':img.select('red'),
          'green':img.select('green')
        });
  
      var ndsi = img.expression(
        '((swir+red)-(blue+nir))/(swir+red+blue+nir)',
        {
          'swir':img.select('swir1'),
          'nir':img.select('nir1'),
          'red':img.select('red'),
          'green':img.select('green'),
          'blue':img.select('blue')
        });
        
      var ndbsi=ee.Image(ibi).add(ee.Image(ndsi)).divide(2)
      return img.addBands(ndsi.rename("NDBSI"));
    }
}

/************************************************************************
 * 
 *                            normalization
 * 
 * *********************************************************************/



//============================================ step1 filter water
var al_single=col_a1.median().select(['sur_refl_b01','sur_refl_b02','sur_refl_b03','sur_refl_b04','sur_refl_b05','sur_refl_b06','sur_refl_b07'],
                                     ["red","nir1","blue","green","nir2","swir1","swir2"]).clip(roi);
var mndwi=a1_index.mNDWI(al_single).select("mNDWI")
var water_mask = mndwi.lt(0.2);
Map.addLayer(water_mask,{palette:['ffffff','000000']},'water mask');

//============================================ step2 mask water
al_single=al_single.updateMask(water_mask);

//============================================ step3 calculate index
var lst=col_lst.mean().updateMask(water_mask).clip(roi).rename('LST')
var ndvi=a1_index.NDVI(al_single).select("NDVI")
var wet=a1_index.WET(al_single).select("WET")
var ndbsi=a1_index.NDBSI(al_single).select("NDBSI")

//============================================ step4 combine
var bands_4=ee.Image.cat([ndvi,wet,ndbsi,lst]) 


Map.addLayer(wet,{min:-0.5,max:0.6,palette:['adff2f','c0e366', 'ccc68b', 'd2a8aa', 'd388c7', 'd063e3', 'c72fff']},'wet');
Map.addLayer(lst,{min:10,max:40,palette:['adff2f','c0e366', 'ccc68b', 'd2a8aa', 'd388c7', 'd063e3', 'c72fff']},'LST');
Map.addLayer(ndvi,{min:-1,max:1,palette:['ff4b2f','ffc0cb', 'd388c7', 'adff2f', '7dd100']},'ndvi');
Map.addLayer(ndbsi,{min:-1,max:1,palette:['adff2f','c0e366', 'ccc68b', 'd2a8aa', 'd388c7', 'd063e3', 'c72fff']},'ndbsi');


//============================================ step5 normalization


//标准化 
function normalization(image,region,scale){
   
   var num = image.reduceRegion({
        reducer:ee.Reducer.percentile([1,99]),
        geometry:region,
        scale:scale,
        maxPixels:1e13,
        // tileScale: 16
      })
  
// use unit scale to normalize the pixel values
  var unitScale = ee.ImageCollection.fromImages(
    image.bandNames().map(function(name){
    name = ee.String(name);
    var num_1 = ee.Number(num.get(name.cat('_p1'))); // get the minimum cutoff value
    var num_99 = ee.Number(num.get(name.cat('_p99'))); // get the maximum cutoff value
    var band = image.select(name);
    var max = num_99;
    var min = num_1;
    var band1=ee.Image(min).multiply(band.lt(min)).add(ee.Image(max).multiply(band.gt(max)))
                        .add(band.multiply(ee.Image(1).subtract(band.lt(min)).subtract(band.gt(max))))
    var result_band=band1.subtract(min).divide(max.subtract(min));
    return result_band;
  })).toBands().rename(image.bandNames()); // conver to band
    return unitScale;
}

var normal_img=normalization(bands_4,roi,scale_my)


//============================================ step6 calculate mean
// calculate the mean by origain image 
var nor_mal = normal_img.select(['NDBSI','WET','NDVI','LST']).reduceRegion({
    reducer: ee.Reducer.mean()
              .combine(ee.Reducer.stdDev(),null, true),
    geometry: roi,
    scale: scale_my,
    maxPixels: 1e13,
    // tileScale: 16
  }); 

// calculate the mean by normalization image 
 var roi_mal = bands_4.select(['NDBSI','WET','NDVI','LST']).reduceRegion({
    reducer: ee.Reducer.mean()
              .combine(ee.Reducer.stdDev(),null, true),
    geometry: roi,
    scale: scale_my,
    maxPixels: 1e13,
    // tileScale: 16
  }); 

print('ori_mean',roi_mal)

print('normal_mean',nor_mal)


/************************************************************************
 * 
 *                            PCA
 * 
 * *********************************************************************/
 print('============================================> PCA ')
// get pca
var flag=ee.Array([1,1,-1,-1])
function pca_model(image,scale){
  var bandNames = image.bandNames();
  print('bandNames',bandNames)
  var region = roi;
  var meanDict = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry:region,
    scale: scale,
    maxPixels: 1e13});
  var means = ee.Image.constant(meanDict.values(bandNames));
  // print(means);
  var centered = image.subtract(means); // get the de-mean value
  var getNewBandNames = function(prefix) {
  var seq = ee.List.sequence(1, bandNames.length());
  return seq.map(function(b) {
    return ee.String(prefix).cat(ee.Number(b).int());
  })};
  var arrays = centered.toArray(); // centerise and image to array
  var covar = arrays.reduceRegion({
    reducer: ee.Reducer.centeredCovariance(), // get the cov matrix
    geometry: region,
    scale: scale,
    maxPixels: 1e13
  });
  // print(covar);
  var covarArray = ee.Array(covar.get('array'));
  var eigens = covarArray.eigen(); // get the vector
  // print('eigens',eigens);
  // var eigens_list = ee.List(eigens); // convert to the array
  
  // print(eigens_array.transpose(),'eigens_array');
  var eigens_img = ee.Image(eigens).rename('eigien');
  
  var eigenValues = eigens.slice(1, 0, 1);
  // print(eigenValues,'eigen value');
  var eigenVectors = eigens.slice(1, 1);
  // print(eigenVectors,'eigen vector');
  var arrayImage = arrays.toArray(1); // to single image, contain 4 values
  var eigenVectors=ee.Array([ee.Array(eigenVectors.toList().get(0)).abs().multiply(flag).toList(),
                  eigenVectors.toList().get(1),
                  eigenVectors.toList().get(2),
                  eigenVectors.toList().get(3)])
  
  // print(arrayImage);
  // Map.addLayer(arrayImage);
  var principalComponents = ee.Image(eigenVectors).matrixMultiply(arrayImage); // get the final results
  var sdImage = ee.Image(eigenValues.sqrt())
    .arrayProject([0]).arrayFlatten([getNewBandNames('sd')]);
  // print(sdImage);
  // Map.addLayer(sdImage);

 // eigenvalue
  var eigenValues1 = ee.List(eigenValues.toList().get(0)).get(0);
  var eigenValues2 = ee.List(eigenValues.toList().get(1)).get(0);
  var eigenValues3 = ee.List(eigenValues.toList().get(2)).get(0);
  var eigenValues4 = ee.List(eigenValues.toList().get(3)).get(0);
  // var eigenValues5 = ee.List(eigenValues.toList().get(4)).get(0);
  
  var eigen_sum=eigenValues.toList().flatten().reduce(ee.Reducer.sum())
  var contribution_rate=eigenValues.divide(eigen_sum)
  
  print('covarArray',covarArray)
  print('eigenVectors',eigenVectors)
  print('eigenValues',eigenValues)
  print('contribution_rate',contribution_rate)
  
  // eigenvector
  var pca1 = principalComponents.arrayProject([0]).arrayFlatten([getNewBandNames('pc')]).divide(sdImage).select(0);
  var pca2 = principalComponents.arrayProject([0]).arrayFlatten([getNewBandNames('pc')]).divide(sdImage).select(1);
  var pca3 = principalComponents.arrayProject([0]).arrayFlatten([getNewBandNames('pc')]).divide(sdImage).select(2);
  var pca4 = principalComponents.arrayProject([0]).arrayFlatten([getNewBandNames('pc')]).divide(sdImage).select(3);
  
  
  return pca1

}
var pca1 = pca_model(normal_img,scale_my);
//pca1=ee.Image(1).clip(roi).subtract(pca1)
/****************************************************************************************************************************
 * 
 *                                             RSEI
 * 
 * **************************************************************************************************************************/

var rsei = normalization(pca1,roi,scale_my).rename('rsei');
// var raw_chart=ui.Chart.image.histogram(rsei.select('rsei'),xlglm_g,30);
// print(raw_chart);

Map.addLayer(rsei,{palette:['FF4500','FFA500', 'EEE8AA', '228B22','006400'],min:0,max:1},'rsei')
var mean = rsei.reduceRegion({
  reducer:ee.Reducer.mean(),
  geometry:roi,
  scale:scale_my,
  maxPixels: 1e13
})
print(mean,'mean');

/****************************************************************************************************************************
 * 
 *                                             Export
 * 
 * **************************************************************************************************************************/
 var test=normal_img
 Export.image.toDrive({
  image:bands_4.clip(roi),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_bands4_index_original',
  description:year+'_bands4_index_original',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})

Export.image.toDrive({
  image:test.clip(roi),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_bands4_index_normal',
  description:year+'_bands4_index_normal',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})


Export.image.toDrive({
  image:rsei.clip(roi),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_rsei',
  description:year+'_rsei',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})


Export.image.toDrive({
  image:al_single.select(["red","green","blue"]).clip(roi),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_bands_rgb',
  description:year+'_bands_rgb',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})
///////////////////////////////////////////////////////////////
// test.select(['ndvi','wet','LST','ndbsi'])
 Export.image.toDrive({
  image:test.select(['ndvi']).clip(roi),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_ndvi',
  description:year+'_ndvi',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})

// test.select(['ndvi','wet','LST','ndbsi'])
 Export.image.toDrive({
  image:test.select(['wet']).clip(roi),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_wet',
  description:year+'_wet',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})

// test.select(['ndvi','wet','LST','ndbsi'])
 Export.image.toDrive({
  image:test.select(['LST']).clip(roi).subtract(273.15),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_LST',
  description:year+'_LST',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})

// test.select(['ndvi','wet','LST','ndbsi'])
 Export.image.toDrive({
  image:test.select(['ndbsi']).clip(roi),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_ndbsi',
  description:year+'_ndbsi',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})

 Export.image.toDrive({
  image:mndwi.gte(0.2),
  folder:'Rseis_modis',
  fileNamePrefix :year+'_water',
  description:year+'_water',
  region:roi,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
})
