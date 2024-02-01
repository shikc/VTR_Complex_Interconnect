#ifndef RR_GRAPH2_H
#define RR_GRAPH2_H
#include <vector>

#include "build_switchblocks.h"
#include "rr_graph_util.h"
#include "rr_types.h"
#include "device_grid.h"

/******************* Types shared by rr_graph2 functions *********************/

enum e_seg_details_type {
    SEG_DETAILS_X,
    SEG_DETAILS_Y
};

struct t_rr_edge_info {
    t_rr_edge_info(int from, int to, short type) noexcept
        : from_node(from)
        , to_node(to)
        , switch_type(type) {}

    int from_node = OPEN;
    int to_node = OPEN;
    short switch_type = OPEN;

    friend bool operator<(const t_rr_edge_info& lhs, const t_rr_edge_info& rhs) {
        return std::tie(lhs.from_node, lhs.to_node, lhs.switch_type) < std::tie(rhs.from_node, rhs.to_node, rhs.switch_type);
    }

    friend bool operator==(const t_rr_edge_info& lhs, const t_rr_edge_info& rhs) {
        return std::tie(lhs.from_node, lhs.to_node, lhs.switch_type) == std::tie(rhs.from_node, rhs.to_node, rhs.switch_type);
    }
};

typedef std::vector<t_rr_edge_info> t_rr_edge_info_set;

typedef vtr::NdMatrix<short, 6> t_sblock_pattern;

/******************* Subroutines exported by rr_graph2.c *********************/
t_medium_rr_node_indices alloc_and_load_medium_rr_node_indices(const DeviceGrid& grid,
                                                                int* index,
                                                                const t_arch* arch,
                                                                const t_medium_seg_inf& medium_seg_inf,
                                                                const t_chan_details& chan_details_x,
                                                                const t_chan_details& chan_details_y);

t_rr_node_indices alloc_and_load_rr_node_indices(const int max_chan_width,
                                                 const DeviceGrid& grid,
                                                 int* index,
                                                 const t_chan_details& chan_details_x,
                                                 const t_chan_details& chan_details_y);

//Returns all the rr nodes associated with the specified coordinate (i.e. accross sides)
std::vector<int> get_rr_node_indices(const t_rr_node_indices& L_rr_node_indices,
                                     int x,
                                     int y,
                                     t_rr_type rr_type,
                                     int ptc);

//Returns all x-channel or y-channel wires at the specified location
std::vector<int> get_rr_node_chan_wires_at_location(const t_rr_node_indices& L_rr_node_indices,
                                                    t_rr_type rr_type,
                                                    int x,
                                                    int y);

int find_average_rr_node_index(int device_width,
                               int device_height,
                               t_rr_type rr_type,
                               int ptc,
                               const t_rr_node_indices& L_rr_node_indices);

t_seg_details* alloc_and_load_seg_details(int* max_chan_width,
                                          const int max_len,
                                          const t_arch* arch,
                                          t_medium_seg_inf& medium_seg_inf,
                                          const std::vector<t_segment_inf>& segment_inf,
                                          const bool use_full_seg_groups,
                                          const bool is_global_graph,
                                          const enum e_directionality directionality,
                                          std::vector<int>& nums_per_segtype,
                                          std::map<int, int>& clbpin_medium_map,
                                          std::vector<std::map<std::string, int>>& mux_name_medium_map_vec,
                                          bool& use_bend_segment_groups,
                                          int* num_seg_details = nullptr);

void alloc_and_load_chan_details(const DeviceGrid& grid,
                                 const t_chan_width* nodes_per_chan,
                                 const bool trim_empty_channels,
                                 const bool trim_obs_channels,
                                 const int num_seg_details,
                                 const t_seg_details* seg_details,
                                 const t_medium_seg_inf medium_seg_inf,
                                 t_chan_details& chan_details_x,
                                 t_chan_details& chan_details_y);
t_chan_details init_chan_details(const DeviceGrid& grid,
                                 const t_chan_width* nodes_per_chan,
                                 const int num_seg_details,
                                 const t_seg_details* seg_details,
                                 const t_medium_seg_inf medium_seg_inf,
                                 const enum e_seg_details_type seg_details_type);
void obstruct_chan_details(const DeviceGrid& grid,
                           const t_chan_width* nodes_per_chan,
                           const t_medium_seg_inf medium_seg_inf,
                           const bool trim_empty_channels,
                           const bool trim_obs_channels,
                           t_chan_details& chan_details_x,
                           t_chan_details& chan_details_y);
void adjust_chan_details(const DeviceGrid& grid,
                         const t_chan_width* nodes_per_chan,
                         t_chan_details& chan_details_x,
                         t_chan_details& chan_details_y);
void adjust_seg_details(const int x,
                        const int y,
                        const DeviceGrid& grid,
                        const t_chan_width* nodes_per_chan,
                        t_chan_details& chan_details,
                        const enum e_seg_details_type seg_details_type);

void free_chan_details(t_chan_details& chan_details_x,
                       t_chan_details& chan_details_y);

int get_seg_start(const t_chan_seg_details* seg_details,
                  const int itrack,
                  const int chan_num,
                  const int seg_num);
int get_seg_end(const t_chan_seg_details* seg_details,
                const int itrack,
                const int istart,
                const int chan_num,
                const int seg_max);

bool is_cblock(const int chan,
               const int seg,
               const int track,
               const t_chan_seg_details* seg_details);

bool is_sblock(const int chan,
               int wire_seg,
               const int sb_seg,
               const int track,
               const t_chan_seg_details* seg_details,
               const enum e_directionality directionality);

int get_bidir_opin_connections(const int i,
                               const int j,
                               const int ipin,
                               const int from_rr_node,
                               t_rr_edge_info_set& rr_edges_to_create,
                               const t_pin_to_track_lookup_bi& opin_to_track_map,
                               const t_rr_node_indices& L_rr_node_indices,
                               const t_chan_details& chan_details_x,
                               const t_chan_details& chan_details_y);

int get_unidir_opin_connections_inc(const int chan,
                                    const int seg,
                                    int Fc,
                                    const int seg_type_index,
                                    const t_rr_type chan_type,
                                    const t_chan_seg_details* seg_details,
                                    const int from_rr_node,
                                    t_rr_edge_info_set& rr_edges_to_create,
                                    vtr::NdMatrix<int, 4>& Fc_ofs,
                                    const int max_len,
                                    const int max_chan_width,
                                    const t_rr_node_indices& L_rr_node_indices,
                                    bool* Fc_clipped,
                                    const DeviceGrid& grid);

int get_unidir_opin_connections_dec(const int chan,
                                    const int seg,
                                    int Fc,
                                    const int seg_type_index,
                                    const t_rr_type chan_type,
                                    const t_chan_seg_details* seg_details,
                                    const int from_rr_node,
                                    t_rr_edge_info_set& rr_edges_to_create,
                                    vtr::NdMatrix<int, 4>& Fc_ofs,
                                    const int max_len,
                                    const int max_chan_width,
                                    const t_rr_node_indices& L_rr_node_indices,
                                    bool* Fc_clipped,
                                    const DeviceGrid& grid);

int get_track_to_pins(int seg,
                      int chan,
                      int track,
                      int tracks_per_chan,
                      int from_rr_node,
                      t_rr_edge_info_set& rr_edges_to_create,
                      const t_rr_node_indices& L_rr_node_indices,
                      const t_track_to_pin_lookup& track_to_pin_lookup,
                      const t_chan_seg_details* seg_details,
                      enum e_rr_type chan_type,
                      int chan_length,
                      int wire_to_ipin_switch,
                      enum e_directionality directionality);

int get_track_to_tracks(const int from_chan,
                        const int from_seg,
                        const int from_track,
                        const t_rr_type from_type,
                        const int to_seg,
                        const t_rr_type to_type,
                        const int chan_len,
                        const int max_chan_width,
                        const DeviceGrid& grid,
                        const int Fs_per_side,
                        t_sblock_pattern& sblock_pattern,
                        const int from_rr_node,
                        t_rr_edge_info_set& rr_edges_to_create,
                        const t_chan_seg_details* from_seg_details,
                        const t_chan_seg_details* to_seg_details,
                        const t_chan_details& to_chan_details,
                        const enum e_directionality directionality,
                        const t_rr_node_indices& L_rr_node_indices,
                        const vtr::NdMatrix<std::vector<int>, 3>& switch_block_conn,
                        t_sb_connection_map* sb_conn_map,
                        const int bend_delayless_switch,
                        std::map<int, int>& bend_segment_map);

int get_track_to_tracks_normal_when_gsb(const int from_chan,
                                        const int from_seg,
                                        const int from_track,
                                        const t_rr_type from_type,
                                        const int to_seg,
                                        const t_rr_type to_type,
                                        const int sb_seg,
                                        const int chan_len,
                                        const int max_chan_width,
                                        const DeviceGrid& grid,
                                        const int Fs_per_side,
                                        t_sblock_pattern& sblock_pattern,
                                        const int from_rr_node,
                                        t_rr_edge_info_set& rr_edges_to_create,
                                        const t_chan_seg_details* from_seg_details,
                                        const t_chan_seg_details* to_seg_details,
                                        const t_chan_details& to_chan_details,
                                        const enum e_directionality directionality,
                                        const t_rr_node_indices& L_rr_node_indices,
                                        const int bend_delayless_switch,
                                        std::map<int, int>& bend_segment_map);

t_sblock_pattern alloc_sblock_pattern_lookup(const DeviceGrid& grid,
                                             const int max_chan_width);
                                             
void load_sblock_bend_pattern_lookup_for_gsb(const int i,
                                             const int j,
                                             const DeviceGrid& grid,
                                             const t_chan_width* nodes_per_chan,
                                             const t_chan_details& chan_details_x,
                                             const t_chan_details& chan_details_y,
                                             const int Fs,
                                             const enum e_switch_block_type switch_block_type,
                                             t_sblock_pattern& sblock_pattern);

void load_sblock_pattern_lookup(const int i,
                                const int j,
                                const DeviceGrid& grid,
                                const t_chan_width* nodes_per_chan,
                                const t_chan_details& chan_details_x,
                                const t_chan_details& chan_details_y,
                                const int Fs,
                                const enum e_switch_block_type switch_block_type,
                                t_sblock_pattern& sblock_pattern);

int* get_seg_track_counts(const int num_sets,
                          const std::vector<t_segment_inf>& segment_inf,
                          const bool use_full_seg_groups);

void dump_seg_details(const t_chan_seg_details* seg_details,
                      int max_chan_width,
                      const char* fname);
void dump_seg_details(const t_chan_seg_details* seg_details,
                      int max_chan_width,
                      FILE* fp);
void dump_chan_details(const t_chan_details& chan_details_x,
                       const t_chan_details& chan_details_y,
                       int max_chan_width,
                       const DeviceGrid& grid,
                       const char* fname);
void dump_sblock_pattern(const t_sblock_pattern& sblock_pattern,
                         int max_chan_width,
                         const DeviceGrid& grid,
                         const char* fname);

void dump_custom_sblock_pattern(t_sb_connection_map* sb_conn_map,
                                int max_chan_width,
                                const DeviceGrid& grid,
                                const char* fname);

//Partitions RR graph edges to allow fast access to configurable/non-configurabe edge subsets
void partition_rr_graph_edges(DeviceContext& device_ctx);


//for gsb arch
int get_track_to_track_for_gsb(int x_coord, int y_coord, const DeviceGrid& grid,
                                        const t_chan_details& chan_details_x,
                                        const t_chan_details& chan_details_y,
                                        const t_rr_node_indices& L_rr_node_indices,
                                        t_rr_edge_info_set& rr_edges_to_create,
                                        t_gsb_connection_map* gsb_conn_map,
                                        t_gsb_multistage_mux_map* gsb_multistage_mux_map);

int get_omuxOrpb_to_track_for_gsb(int x_coord, int y_coord, const DeviceGrid& grid,
                                    const t_chan_details& chan_details_x,
                                    const t_chan_details& chan_details_y,
                                    const t_rr_node_indices& L_rr_node_indices,
                                    const t_medium_rr_node_indices& L_medium_rr_node_indices,
                                    t_rr_edge_info_set& rr_edges_to_create,
                                    t_gsb_connection_map* gsb_conn_map,
                                    const std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map,
                                    t_gsb_multistage_mux_map* gsb_multistage_mux_map,
                                    e_from_type from_type);

int build_multistage_mux_for_gsb(int x_coord, int y_coord, t_rr_edge_info_set& rr_edges_to_create, t_gsb_multistage_mux_map* gsb_multistage_mux_map);

int get_track_to_imux_for_imux(int x_coord, int y_coord, const DeviceGrid& grid,
                                    const t_chan_details& chan_details_x,
                                    const t_chan_details& chan_details_y,
                                    const t_rr_node_indices& L_rr_node_indices,
                                    const t_medium_rr_node_indices& L_medium_rr_node_indices,
                                    t_rr_edge_info_set& rr_edges_to_create,
                                    t_gsb_connection_map* imux_conn_map,
                                    const std::map<int, int>& clbpin_to_medium,
                                    const std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map);

int get_omuxOrpb_to_imux_for_imux(int x_coord, int y_coord,
                                    const std::map<std::string, std::pair<int, int>>& record_imux_coord, 
                                    const DeviceGrid& grid,
                                    const t_rr_node_indices& L_rr_node_indices,
                                    const t_medium_rr_node_indices& L_medium_rr_node_indices,
                                    t_rr_edge_info_set& rr_edges_to_create,
                                    t_gsb_connection_map* imux_conn_map,
                                    const std::map<int, int>& clbpin_to_medium,
                                    const std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map,
                                    std::string actual_pbtype_name,
                                    e_from_type from_type);

int get_medium_to_pin_for_imux(int x_coord, int y_coord,
                                  const std::map<std::string, std::pair<int, int>>& record_imux_coord,
                                  const DeviceGrid& grid,
                                  const t_rr_node_indices& L_rr_node_indices,
                                  const t_medium_rr_node_indices& L_medium_rr_node_indices,
                                  t_rr_edge_info_set& rr_edges_to_create,
                                  const std::map<int, int>& clbpin_to_medium,
                                  const std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map,
                                  std::string actual_pbtype_name,
                                  const int delayless_switch);

int get_pb_to_omux(int x_coord, int y_coord,
                    const DeviceGrid& grid,
                    const t_rr_node_indices& L_rr_node_indices,
                    const t_medium_rr_node_indices& L_medium_rr_node_indices,
                    t_rr_edge_info_set& rr_edges_to_create,
                    t_omux_connection_map* omux_conn_map,
                    std::string actual_pbtype_name,
                    const std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map);

void build_other_pbpin_to_clbpin(const t_gsb_inf& gsb, std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map);

int find_average_medium_rr_node_index(int device_width,
                                      int device_height,
                                      int medium_type,
                                      int ptc,
                                      const t_medium_rr_node_indices& L_medium_rr_node_indices);

int connect_bounds_wires(const t_wire_type_sizes& wire_type_sizes, 
                                const DeviceGrid& grid,
                                const t_chan_details& chan_details_x,
                                const t_chan_details& chan_details_y,
                                const t_rr_node_indices& L_rr_node_indices,
                                t_rr_edge_info_set& rr_edges_to_create,
                                const int delayless_switch);

int build_two_stage_mux_for_gsb(int x_coord, int y_coord, t_rr_edge_info_set& rr_edges_to_create, t_two_stage_mux_map* two_stage_mux_map);

void filter_midpoint_edges(t_rr_edge_info_set& rr_edges_to_create, const std::vector<t_segment_inf>& segment_inf);
#endif
