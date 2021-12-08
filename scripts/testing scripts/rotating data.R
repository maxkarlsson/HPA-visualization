

apa <- 
  get_edges()(g$data) %>%
  as_tibble() %>%
  left_join(color_mapping %>%
              select(label = label_col,
                     color = color_col),
            by = c("node2.label" = "label"))

x <- apa$x
y <- apa$y

tibble(x, y) %>% 
  ggplot(aes(x, y)) +
  geom_point() +
  geom_path()


rotate_coords <- 
  function(x, y, rotate_angle, rotate_center = c(0, 0)) {
    
    # Center data
    rotdata <- 
      tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             
             # Angle
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)) + rotate_angle,
             
             # New coordinates
             x = cos(angle) * hyp,
             y = sin(angle) * hyp,
             
             # Recenter coordinates
             x = x - rotate_center[1],
             y = y - rotate_center[2])
    
    
    rotdata
  }
  
rot_data <- 
  tibble(rotate_angle = seq(0, pi, length.out = 100)) %>% 
  group_by(rotate_angle) %>% 
  do({
    rotate_coords(x, y, .$rotate_angle)
  })

rot_data_offcenter <- 
  tibble(rotate_angle = seq(0, 3, length.out = 10)) %>% 
  group_by(rotate_angle) %>% 
  do({
    rotate_coords(x, y, .$rotate_angle, rotate_center = c(0.5, 0))
  })


library(gganimate)
plot <- 
  rot_data %>% 
  mutate(rotate_angle = factor(rotate_angle, unique(rotate_angle))) %>% 
  ggplot(aes(x, y)) +
  geom_point() +
  # geom_path() +
  transition_states(rotate_angle,
                    transition_length = 2,
                    state_length = 0)


anim <- 
  animate(plot,
          width = 300, height = 300, res = 50,
          nframes = 50)
anim


anim_save(savepath("rotating_points.gif"), anim)

plot <- 
  rot_data_offcenter %>% 
  mutate(rotate_angle = factor(rotate_angle, unique(rotate_angle))) %>% 
  ggplot(aes(x, y)) +
  geom_point() +
  # geom_path() +
  transition_states(rotate_angle,
                    transition_length = 2,
                    state_length = 0)


anim <- 
  animate(plot,
          width = 300, height = 300, res = 50,
          nframes = 50)
anim

######


shrink_rotation_coords <- 
  function(x, y, shrink_angle, rotate_center = c(0, 0)) {
    
    
    shrink_factor <- 
      1 - shrink_angle / (2 * pi)
    
    # Center data
    rotdata <- 
      tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             
             # Angle
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp)),
             
             # Shrink angle
             # angle = case_when(x == 0 & y == 0 ~ 0, 
             #                   quadrant %in% c(1,4) ~ angle + shrink_angle,
             #                   quadrant %in% c(2,3) ~ angle - shrink_angle),
             
             angle = angle * shrink_factor,
             
             # New coordinates
             x = cos(angle) * hyp,
             y = sin(angle) * hyp,
             
             # Recenter coordinates
             x = x - rotate_center[1],
             y = y - rotate_center[2])
    
    
    rotdata 
  }


calculate_coord_angle <- 
  function(x, y, rotate_center = c(0, 0)) {
    tibble(x, y) %>% 
      mutate(x = x + rotate_center[1],
             y = y + rotate_center[2],
             
             
             
             # Calculate quadrants
             quadrant = case_when(x >= 0 & y >= 0 ~ 1,
                                  x < 0 & y >= 0 ~ 2,
                                  x < 0 & y < 0 ~ 3,
                                  x >= 0 & y < 0 ~ 4),
             
             
             
             # Hypotenuse
             hyp = sqrt(x^2 + y^2),
             
             
             # Angle
             
             angle = case_when(x == 0 & y == 0 ~ 0, 
                               quadrant %in% 1:2 ~ acos(x/hyp),
                               quadrant %in% 3:4 ~ 2 * pi - acos(x/hyp))) %>% 
      pull(angle)
  }


calculate_retina_cut_angle <- 
  function(clust) {
    dendrogram <-
      clust %>%
      as.dendrogram()
    
    g <-
      ggraph(dendrogram, layout = 'dendrogram', circular = T)
    
    g_edgepoints <- 
      g$data %>% 
      as_tibble() %>% 
      filter(height == 0) %>% 
      left_join(cutree(clust, k = 2) %>% 
                  enframe("label", 
                          "cluster"),
                by = "label") %>% 
      mutate(angle = calculate_coord_angle(x, y))
    
    expand_grid(node1 = g_edgepoints$.ggraph.index,
                node2 = g_edgepoints$.ggraph.index) %>% 
      left_join(g_edgepoints %>% 
                  select(node1 = .ggraph.index,
                         angle1 = angle, 
                         cluster1 = cluster),
                by = "node1") %>% 
      left_join(g_edgepoints %>% 
                  select(node2 = .ggraph.index,
                         angle2 = angle, 
                         cluster2 = cluster),
                by = "node2") %>% 
      filter(cluster1 == 1,
             cluster2 == 2) %>% 
      group_by_all() %>% 
      mutate(dist = c(angle1 - angle2,
                      (angle1 - 2 * pi) - angle2,
                      angle1 - (angle2 - 2 * pi),
                      (angle1 - 2 * pi) - (angle2 - 2 * pi)) %>% 
               abs() %>% 
               min()) %>% 
      ungroup() %>% 
      arrange(dist) %>% 
      slice(1:2) %>% 
      mutate(cut_angle = (angle1 + angle2) / 2) %>% 
      pull(cut_angle)
  }
