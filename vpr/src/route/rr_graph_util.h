#ifndef RR_GRAPH_UTIL_H
#define RR_GRAPH_UTIL_H

int seg_index_of_cblock(t_rr_type from_rr_type, int to_node);

int seg_index_of_sblock(int from_node, int to_node);

//Return the first rr node of the specified type and coordinates
// For non-IPIN/OPIN types 'side' is ignored
int get_rr_node_index(const t_rr_node_indices& L_rr_node_indices,
                      int x,
                      int y,
                      t_rr_type rr_type,
                      int ptc,
                      e_side side = NUM_SIDES);

int get_medium_rr_node_index(const t_medium_rr_node_indices& L_medium_rr_node_indices,
                             int x,
                             int y,
                             int medium_type,  //0-IMUX, 1-OMUX, 2-GSB
                             int num_type_ind, //类型索引，如今默认都只有一种
                             int ptc);

#endif
