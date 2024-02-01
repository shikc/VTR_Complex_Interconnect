#ifndef RR_GRAPH_INDEXED_DATA_H
#define RR_GRAPH_INDEXED_DATA_H
#include "physical_types.h"

void alloc_and_load_rr_indexed_data(const std::vector<t_segment_inf>& segment_inf,
                                    const t_medium_seg_inf& medium_seg_inf,
                                    const t_rr_node_indices& L_rr_node_indices,
                                    const t_medium_rr_node_indices& L_medium_rr_node_indices,
                                    const int nodes_per_chan,
                                    int wire_to_ipin_switch,
                                    enum e_base_cost_type base_cost_type);

void load_rr_index_segments(const int num_segment, const std::vector<t_segment_inf>& segment_inf);

#endif
