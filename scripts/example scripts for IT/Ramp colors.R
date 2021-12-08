


sc_colors <- 
  read_tsv("../../Data/HPA/HPA21_E103_files/colors/Single cell type.txt",col_names = c("group", "cell", "color"))

biases <- seq(0, 0.8, length.out = 21)
sc_colors %>% 
  select(group, color) %>% 
  distinct() %>% 
  group_by(group) %>% 
  summarise(color = ramp_color_bias(color, "white", biases)) %>% 
  mutate(i = row_number()) %>% 
  ungroup() %>% 
  ggplot(aes(i, group, fill = color)) +
  geom_tile() +
  scale_fill_identity() +
  theme_void() +
  theme(axis.text.y = element_text(hjust = 1))
