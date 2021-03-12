queries =  { 'bipartite_tumor_stroma' : 
 
 '''SELECT
  distance,
  cell_id_1,
  cell_id_2,
  tissue_category_1,
  tissue_category_2
FROM (
  SELECT
    table_1.cell_id_1,
    table_2.cell_id_2,
    table_1.x_1,
    table_1.y_1,
    table_2.x_2,
    table_2.y_2,
    tissue_category_1,
    tissue_category_2 ,

    SQRT(POW((x_1-x_2),2)+POWER((y_1-y_2),2)) AS distance
  FROM (
    SELECT
      cell_x_position AS x_1,
      cell_y_position AS y_1,
      cell_id AS cell_id_1,
      tissue_category as tissue_category_1
    FROM
      `advance-sonar-306410.deepmelo.cell_data`
      Where (tissue_category = 'stroma' )) AS table_1
  CROSS JOIN (
    SELECT
      cell_x_position AS x_2,
      cell_y_position AS y_2,
      cell_id AS cell_id_2,
      tissue_category as tissue_category_2
    FROM
      `advance-sonar-306410.deepmelo.cell_data` 
      Where (tissue_category = 'tumor')) AS table_2 )''',
 
 
 #####################################
 
 
 'distance_statistics_all_cells' : 
 
 '''SELECT 
 
cell_id_1, 
AVG(distance) as mean_distance,
VARIANCE(distance) as var_distance,
MAX(distance) as max_distance,
min(distance) as min_distance




FROM (
SELECT
  distance,
  cell_id_1,
  cell_id_2
  tissue_category_1,
  tissue_category_2
FROM (
  SELECT
    table_1.cell_id_1,
    table_2.cell_id_2,
    table_1.x_1,
    table_1.y_1,
    table_2.x_2,
    table_2.y_2,
    tissue_category_1,
    tissue_category_2 ,
    SQRT(POW((x_1-x_2),2)+POWER((y_1-y_2),2)) AS distance

  FROM (
    SELECT
      cell_x_position AS x_1,
      cell_y_position AS y_1,
      cell_id AS cell_id_1,
      tissue_category as tissue_category_1
    FROM
      `advance-sonar-306410.deepmelo.cell_data`
      Where (cell_id>=0)  ) AS table_1
  CROSS JOIN (
    SELECT
      cell_x_position AS x_2,
      cell_y_position AS y_2,
      cell_id AS cell_id_2,
      tissue_category as tissue_category_2
    FROM
      `advance-sonar-306410.deepmelo.cell_data` 
      Where (cell_id>0)) AS table_2 )

   WHERE cell_id_1 <> cell_id_2

       )


GROUP BY  ( cell_id_1) ''' ,
 
 
 ###############################################################

 
 'distance_statistics_to_stroma' :
 
 
 '''SELECT 
cell_id_1, 
AVG(distance) as mean_distance_to_stroma,
VARIANCE(distance ) as var_distance_to_stroma,
MAX(distance) as max_distance_to_stroma,
min(distance) as min_distance_to_stroma,




FROM (
SELECT
  distance,
  cell_id_1,
  cell_id_2
  tissue_category_1,
  tissue_category_2
FROM (
  SELECT
    table_1.cell_id_1,
    table_2.cell_id_2,
    table_1.x_1,
    table_1.y_1,
    table_2.x_2,
    table_2.y_2,
    tissue_category_1,
    tissue_category_2 ,
    SQRT(POW((x_1-x_2),2)+POWER((y_1-y_2),2)) AS distance

  FROM (
    SELECT
      cell_x_position AS x_1,
      cell_y_position AS y_1,
      cell_id AS cell_id_1,
      tissue_category as tissue_category_1
    FROM
      `advance-sonar-306410.deepmelo.cell_data`
      Where (cell_id>=0)  ) AS table_1
  CROSS JOIN (
    SELECT
      cell_x_position AS x_2,
      cell_y_position AS y_2,
      cell_id AS cell_id_2,
      tissue_category as tissue_category_2
    FROM
      `advance-sonar-306410.deepmelo.cell_data` 
      Where (tissue_category = 'stroma')) AS table_2 )

   WHERE cell_id_1 <> cell_id_2

       )


GROUP BY  ( cell_id_1)   

''' ,
 
####################################################################
 
 
 'distance_statistics_to_tumor' :
 
 
 '''SELECT 
cell_id_1, 
AVG(distance) as mean_distance_to_tumor,
VARIANCE(distance ) as var_distance_to_tumor,
MAX(distance) as max_distance_to_tumor,
min(distance) as min_distance_to_tumor,




FROM (
SELECT
  distance,
  cell_id_1,
  cell_id_2
  tissue_category_1,
  tissue_category_2
FROM (
  SELECT
    table_1.cell_id_1,
    table_2.cell_id_2,
    table_1.x_1,
    table_1.y_1,
    table_2.x_2,
    table_2.y_2,
    tissue_category_1,
    tissue_category_2 ,
    SQRT(POW((x_1-x_2),2)+POWER((y_1-y_2),2)) AS distance

  FROM (
    SELECT
      cell_x_position AS x_1,
      cell_y_position AS y_1,
      cell_id AS cell_id_1,
      tissue_category as tissue_category_1
    FROM
      `advance-sonar-306410.deepmelo.cell_data`
      Where (cell_id>=0)  ) AS table_1
  CROSS JOIN (
    SELECT
      cell_x_position AS x_2,
      cell_y_position AS y_2,
      cell_id AS cell_id_2,
      tissue_category as tissue_category_2
    FROM
      `advance-sonar-306410.deepmelo.cell_data` 
      Where (tissue_category = 'tumor')) AS table_2 )

   WHERE cell_id_1 <> cell_id_2

       )


GROUP BY  ( cell_id_1)   

'''
 
 
##############################################################################
 

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
}