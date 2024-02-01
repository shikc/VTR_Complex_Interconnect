#include <cstdio>
#include <string.h>
#include "limits.h"
using namespace std;

#include "vtr_util.h"
#include "vtr_assert.h"
#include "vtr_log.h"
#include "vtr_memory.h"

#include "vpr_types.h"
#include "vpr_error.h"

#include "globals.h"
#include "rr_graph_util.h"
#include "rr_graph2.h"
#include "rr_graph_sbox.h"
#include "read_xml_arch_file.h"
#include "rr_types.h"

constexpr short UN_SET = -1;

/************************** Subroutines local to this module ****************/

static void get_switch_type(bool is_from_sb,
                            bool is_to_sb,
                            short from_node_switch,
                            short to_node_switch,
                            const int switch_override,
                            short switch_types[2]);

static void load_chan_rr_indices(const int max_chan_width,
                                 const int chan_len,
                                 const int num_chans,
                                 const t_rr_type type,
                                 const t_chan_details& chan_details,
                                 t_rr_node_indices& indices,
                                 const DeviceGrid& grid,
                                 int* index);

static void load_block_rr_indices(const DeviceGrid& grid,
                                  t_rr_node_indices& indices,
                                  int* index);






int count_gsb_medium_num_v2(const t_gsb_inf& gsb, std::map<int, int>& clbpin_medium_map, std::map<std::string, int>& mux_name_medium_map);

void load_imux_medium_rr_indices(const DeviceGrid& grid,
                                  const t_arch* arch,
                                  const t_medium_seg_inf& medium_seg_inf,
                                  t_medium_rr_node_indices& indices,
                                  const t_chan_details& chan_details_x,
                                  const t_chan_details& chan_details_y,
                                  int* index);

void load_omux_medium_rr_indices(const DeviceGrid& grid,
                                  const t_arch* arch,
                                  const t_medium_seg_inf& medium_seg_inf,
                                  t_medium_rr_node_indices& indices,
                                  const t_chan_details& chan_details_x,
                                  const t_chan_details& chan_details_y,
                                  int* index);

void load_gsb_medium_rr_indices(const DeviceGrid& grid,
                                  const t_arch* arch,
                                  const t_medium_seg_inf& medium_seg_inf,
                                  t_medium_rr_node_indices& indices,
                                  const t_chan_details& chan_details_x,
                                  const t_chan_details& chan_details_y,
                                  int* index);

static int get_bidir_track_to_chan_seg(const std::vector<int> conn_tracks,
                                       const t_rr_node_indices& L_rr_node_indices,
                                       const int to_chan,
                                       const int to_seg,
                                       const int to_sb,
                                       const t_rr_type to_type,
                                       const t_chan_seg_details* seg_details,
                                       const bool from_is_sblock,
                                       const int from_switch,
                                       const int switch_override,
                                       const enum e_directionality directionality,
                                       const int from_rr_node,
                                       t_rr_edge_info_set& rr_edges_to_create);

static int get_unidir_track_to_chan_seg(const int from_track,
                                        const int to_chan,
                                        const int to_seg,
                                        const int to_sb,
                                        const t_rr_type to_type,
                                        const int max_chan_width,
                                        const DeviceGrid& grid,
                                        const enum e_side from_side,
                                        const enum e_side to_side,
                                        const int Fs_per_side,
                                        t_sblock_pattern& sblock_pattern,
                                        const int switch_override,
                                        const t_rr_node_indices& L_rr_node_indices,
                                        const t_chan_seg_details* seg_details,
                                        bool* Fs_clipped,
                                        const int from_rr_node,
                                        t_rr_edge_info_set& rr_edges_to_create,
                                        const int bend_delayless_switch,
                                        map<int, int>& bend_segment_map);

static int get_track_to_chan_seg(const int from_track,
                                 const int to_chan,
                                 const int to_seg,
                                 const t_rr_type to_chan_type,
                                 const e_side from_side,
                                 const e_side to_side,
                                 const int swtich_override,
                                 const t_rr_node_indices& L_rr_node_indices,
                                 t_sb_connection_map* sb_conn_map,
                                 const int from_rr_node,
                                 t_rr_edge_info_set& rr_edges_to_create);

static int vpr_to_phy_track(const int itrack,
                            const int chan_num,
                            const int seg_num,
                            const t_chan_seg_details* seg_details,
                            const enum e_directionality directionality);

static bool is_corner(const DeviceGrid& grid, int x, int y);

static bool is_perimeter_top(const DeviceGrid& grid, int x, int y);
static bool is_perimeter_right(const DeviceGrid& grid, int x, int y);
static bool is_perimeter_left(const DeviceGrid& grid, int x, int y);
static bool is_perimeter_bottom(const DeviceGrid& grid, int x, int y);
static bool is_style_perimeter(const DeviceGrid& grid, int x, int y, e_side to_side, e_switch_dir_type switch_type);

static bool label_dangling_side(const e_side side, const int x, const int y, const DeviceGrid& grid);

static int *label_wire_bends(const int chan_num, const int seg_num,
                                const t_chan_seg_details * seg_details, const int seg_type_index, const int max_len,
                                const enum e_direction dir, const int max_chan_width,
                                int *num_bend_wires, bool is_dangling_side);


static int* label_wire_muxes(const int chan_num,
                             const int seg_num,
                             const t_chan_seg_details* seg_details,
                             const int seg_type_index,
                             const int max_len,
                             const enum e_direction dir,
                             const int max_chan_width,
                             const bool check_cb,
                             int* num_wire_muxes,
                             int* num_wire_muxes_cb_restricted,
                             bool is_dangling_side);

static int* label_incoming_bend_wires_ending(const int chan_num,
                                                const int seg_num,
                                                const t_chan_seg_details* seg_details,
                                                const int seg_type_index,
                                                const int max_len,
                                                const enum e_direction dir,
                                                const int max_chan_width,
                                                int* num_incoming_bend,
                                                bool is_dangling_side);

static int* label_incoming_wires(const int chan_num,
                                 const int seg_num,
                                 const int sb_seg,
                                 const t_chan_seg_details* seg_details,
                                 const int max_len,
                                 const enum e_direction dir,
                                 const int max_chan_width,
                                 int* num_incoming_wires,
                                 int* num_ending_wires);

static int find_label_of_track(int* wire_mux_on_track,
                               int num_wire_muxes,
                               int from_track);

void dump_seg_details(t_chan_seg_details* seg_details,
                      int max_chan_width,
                      const char* fname);

void count_bend_wire_type_sizes(const t_chan_seg_details* channel, int bend_start, int bend_end, std::vector<Wire_Info>& real_bend_wire_type);

static int build_a_bound_wire(int x_coord, int y_coord, 
                              const t_wire_type_sizes& wire_type_sizes,
                              const t_chan_details& chan_details,
                              const t_rr_node_indices& L_rr_node_indices,
                              t_rr_edge_info_set& rr_edges_to_create,
                              e_direction to_direction,
                              t_rr_type chan_type,
                              const int delayless_switch);
//Returns how the switch type for the switch block at the specified location should be created
//  grid: The device grid
//  from_chan_coord: The horizontal or vertical channel index (i.e. x-coord for CHANY, y-coord for CHANX)
//  from_seg_coord: The horizontal or vertical location along the channel (i.e. y-coord for CHANY, x-coord for CHANX)
//  from_chan_type: The from channel type
//  to_chan_type: The to channel type
static int should_create_switchblock(const DeviceGrid& grid, int from_chan_coord, int from_seg_coord, t_rr_type from_chan_type, t_rr_type to_chan_type);

static bool should_apply_switch_override(int switch_override);

/******************** Subroutine definitions *******************************/

/* This assigns tracks (individually or pairs) to segment types.
 * It tries to match requested ratio. If use_full_seg_groups is
 * true, then segments are assigned only in multiples of their
 * length. This is primarily used for making a tileable unidir
 * layout. The effect of using this is that the number of tracks
 * requested will not always be met and the result will sometimes
 * be over and sometimes under.
 * The pattern when using use_full_seg_groups is to keep adding
 * one group of the track type that wants the largest number of
 * groups of tracks. Each time a group is assigned, the types
 * demand is reduced by 1 unit. The process stops when we are
 * no longer less than the requested number of tracks. As a final
 * step, if we were closer to target before last more, undo it
 * and end up with a result that uses fewer tracks than given. */
//use_full_seg_groups比例好像是按组来算的，例如3，4倍线各50%，100根，则实际两者组数一样，可能差1.
//普通模式下是按比例直接算的，例如3，4倍线freq均为1,若是单向线，设宽度70，则有35对，一个18对，一个17对
int* get_seg_track_counts(const int num_sets,
                          const std::vector<t_segment_inf>& segment_inf,
                          const bool use_full_seg_groups) {
    int* result;
    double* demand;
    int imax, freq_sum, assigned, size;
    double scale, max, reduce;

    result = (int*)vtr::malloc(sizeof(int) * segment_inf.size());
    demand = (double*)vtr::malloc(sizeof(double) * segment_inf.size());

    /* Scale factor so we can divide by any length
     * and still use integers */
    scale = 1;
    freq_sum = 0;
    for (size_t i = 0; i < segment_inf.size(); ++i) {
        scale *= segment_inf[i].length;
        freq_sum += segment_inf[i].frequency;
    }
    reduce = scale * freq_sum;

    /* Init assignments to 0 and set the demand values */
    for (size_t i = 0; i < segment_inf.size(); ++i) {
        result[i] = 0;
        demand[i] = 0.0;
        demand[i] = scale * num_sets * segment_inf[i].frequency;
        if (use_full_seg_groups || segment_inf[i].isbend) {
            demand[i] /= segment_inf[i].length;
        }
    }

    /* Keep assigning tracks until we use them up */
    assigned = 0;
    size = 0;
    imax = 0;

    // We should assigning tracks to the bend segment first cause they may have low radio
    for (int i_seg = 0; i_seg < segment_inf.size(); ++i_seg) {
        if (segment_inf[i_seg].isbend) {
            while (demand[i_seg] > 0) {
                size = segment_inf[i_seg].length;
                demand[i_seg] -= reduce;
                result[i_seg] += size;
                assigned += size;
            }
        }
    }

    while (assigned < num_sets) {
        /* Find current maximum demand */
        max = 0;
        for (size_t i = 0; i < segment_inf.size(); ++i) {
            if (demand[i] > max) {
                imax = i;
                max = demand[i];
            }
        }

        /* Assign tracks to the type and reduce the types demand */
        size = (use_full_seg_groups ? segment_inf[imax].length : 1);
        demand[imax] -= reduce;
        result[imax] += size;
        assigned += size;
    }

    /* Undo last assignment if we were closer to goal without it */
    if ((assigned - num_sets) > (size / 2)) {
        result[imax] -= size;
    }

    /* Free temps */
    if (demand) {
        vtr::free(demand);
        demand = nullptr;
    }

    /* This must be freed by caller */
    return result;
}

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
                                          int* num_seg_details) {
    /* Allocates and loads the seg_details data structure.  Max_len gives the   *
     * maximum length of a segment (dimension of array).  The code below tries  *
     * to:                                                                      *
     * (1) stagger the start points of segments of the same type evenly;        *
     * (2) spread out the limited number of connection boxes or switch boxes    *
     *     evenly along the length of a segment, starting at the segment ends;  *
     * (3) stagger the connection and switch boxes on different long lines,     *
     *     as they will not be staggered by different segment start points.     */

    int cur_track, ntracks, itrack, length, j, index;
    int arch_wire_switch, arch_opin_switch, fac, num_sets, tmp;
    int group_start, first_track;
    int* sets_per_seg_type = nullptr;
    t_seg_details* seg_details = nullptr;
    bool longline;
    int tmp_sum_len;
    int index_offset;

    bool is_gsb_medium;
    std::string medium_seg_name;

    /* Unidir tracks are assigned in pairs, and bidir tracks individually */
    if (directionality == BI_DIRECTIONAL) {
        fac = 1;
    } else {
        VTR_ASSERT(directionality == UNI_DIRECTIONAL);
        fac = 2;
    }

    if (*max_chan_width % fac != 0) {
        VPR_THROW(VPR_ERROR_ROUTE, "Routing channel width must be divisible by %d (channel width was %d)", fac, *max_chan_width);
    }

    /* Map segment type fractions and groupings to counts of tracks */
    sets_per_seg_type = get_seg_track_counts((*max_chan_width / fac),
                                             segment_inf, use_full_seg_groups);

    /* Count the number tracks actually assigned. */
    tmp = 0;
    for (size_t i = 0; i < segment_inf.size(); ++i) {
        tmp += sets_per_seg_type[i] * fac;
        nums_per_segtype[i] = sets_per_seg_type[i] * fac;
        if (segment_inf[i].isbend)
            use_bend_segment_groups = true;
    }

    VTR_LOG("The segment distribution: \n");
    for (size_t i = 0; i < segment_inf.size(); ++i) {
        VTR_LOG("Segment name : %s Distribution: %d\n", segment_inf[i].name.c_str(), sets_per_seg_type[i] * fac);
    }

    VTR_ASSERT(use_full_seg_groups || (tmp == *max_chan_width));
    *max_chan_width = tmp;

    clbpin_medium_map.clear();
    mux_name_medium_map_vec.resize(1);
    int gsb_medium_num;

    if (arch->gsb_inf[0].twostage_mux_spec) {
        gsb_medium_num = count_gsb_medium_num_v2(arch->gsb_inf[0], clbpin_medium_map, mux_name_medium_map_vec[0]); 
    }

    int medium_num_xory=int(ceil(gsb_medium_num/2.0));
    
    
    std::string medium_seg_names="gsb_medium";
    int arch_switch = segment_inf[segment_inf.size()-1].arch_wire_switch;

    
    medium_seg_inf.gsb_mdeium_num_xory = medium_num_xory;
    medium_seg_inf.gsb_medium_num = gsb_medium_num;

    int total_seg_details = *max_chan_width + medium_num_xory;
    

    seg_details = new t_seg_details[total_seg_details];

    /* Setup the seg_details data */
    cur_track = 0;
    index_offset = 0;
    int cost_index_num = 0;
    for (size_t i = 0; i < segment_inf.size(); ++i) {
        first_track = cur_track;

        num_sets = sets_per_seg_type[i];
        ntracks = fac * num_sets;
        if (ntracks < 1) {
            continue;
        }
        /* Avoid divide by 0 if ntracks */
        longline = segment_inf[i].longline;
        length = segment_inf[i].length;
        if (longline) {
            length = max_len;
        }

        arch_wire_switch = segment_inf[i].arch_wire_switch;
        arch_opin_switch = segment_inf[i].arch_opin_switch;
        VTR_ASSERT((arch_wire_switch == arch_opin_switch) || (directionality != UNI_DIRECTIONAL));

        
        std::vector<int> tmp_seg_len;
        std::vector<e_switch_dir_type> tmp_switch_type;
        if (segment_inf[i].isbend) {
            int tmp_len = 1;
            int sum_len = 0;

            for (int i_len = 0; i_len < segment_inf[i].length - 1; i_len++) {
                if (segment_inf[i].bend[i_len] == 0) {
                    tmp_len++;
                } else if (segment_inf[i].bend[i_len] != 0) {
                    VTR_ASSERT(tmp_len < segment_inf[i].length);
                    tmp_seg_len.push_back(tmp_len);
                    if (segment_inf[i].bend[i_len] == 1)
                        tmp_switch_type.push_back(UP_TYPE);
                    else if (segment_inf[i].bend[i_len] == 2)
                        tmp_switch_type.push_back(DOWN_TYPE);
                    sum_len += tmp_len;
                    tmp_len = 1;
                }
            }
            // add the last clip of segment
            if (sum_len < segment_inf[i].length) {
                tmp_seg_len.push_back(segment_inf[i].length - sum_len);
                tmp_switch_type.push_back(NORMAL_TYPE);
            }
        } else {
            tmp_seg_len.push_back(segment_inf[i].length);
            tmp_switch_type.push_back(NORMAL_TYPE);
        }

        is_gsb_medium = segment_inf[i].is_gsb_medium;

        /* Set up the tracks of same type */
        group_start = 0;
        tmp_sum_len = 0;

        int first_part_wire_switch = -1;
        int first_part_opin_switch = -1;

        for (size_t ipart = 0; ipart < tmp_seg_len.size(); ipart++) {
            for (itrack = 0; itrack < (tmp_seg_len[ipart] * ntracks) / length; itrack++) {
                /* set the name of the segment type this track belongs to */
                seg_details[cur_track].type_name = segment_inf[i].name;
                seg_details[cur_track].part_idx = ipart;
                /* Remember the start track of the current wire group */
                if ((itrack / fac) % tmp_seg_len[ipart] == 0 && (itrack % fac) == 0) {
                    group_start = cur_track;
                }

                seg_details[cur_track].length = tmp_seg_len[ipart];
                seg_details[cur_track].longline = longline;
                seg_details[cur_track].switch_dir_type = tmp_switch_type[ipart];
                seg_details[cur_track].bend_len = tmp_seg_len.size();

                seg_details[cur_track].is_gsb_medium = is_gsb_medium;

                /* Stagger the start points in for each track set. The
             * pin mappings should be aware of this when chosing an
             * intelligent way of connecting pins and tracks.
             * cur_track is used as an offset so that extra tracks
             * from different segment types are hopefully better
             * balanced. */
                seg_details[cur_track].start = (cur_track / fac) % tmp_seg_len[ipart] + 1;

                /* These properties are used for vpr_to_phy_track to determine
             * * twisting of wires. */
                seg_details[cur_track].group_start = group_start;
                seg_details[cur_track].group_size = min(ntracks + first_track - group_start, tmp_seg_len[ipart] * fac);
                VTR_ASSERT(0 == seg_details[cur_track].group_size % fac);
                if (0 == seg_details[cur_track].group_size) {
                    seg_details[cur_track].group_size = tmp_seg_len[ipart] * fac;
                }

                seg_details[cur_track].seg_start = -1;
                seg_details[cur_track].seg_end = -1;

                /* Setup the cb and sb patterns. Global route graphs can't depopulate cb and sb
             * since this is a property of a detailed route. */
                seg_details[cur_track].cb = std::make_unique<bool[]>(tmp_seg_len[ipart]);
                seg_details[cur_track].sb = std::make_unique<bool[]>(tmp_seg_len[ipart] + 1);
                for (j = 0; j < tmp_seg_len[ipart]; ++j) {
                    if (is_global_graph || seg_details[cur_track].longline) {
                        seg_details[cur_track].cb[j] = true;
                    } else {
                        /* Use the segment's pattern. */
                        index = (j + tmp_sum_len) % segment_inf[i].cb.size();
                        seg_details[cur_track].cb[j] = segment_inf[i].cb[index];
                    }
                }
                for (j = 0; j < (tmp_seg_len[ipart] + 1); ++j) {
                    if (is_global_graph || seg_details[cur_track].longline) {
                        seg_details[cur_track].sb[j] = true;
                    } else {
                        /* Use the segment's pattern. */
                        index = (j + tmp_sum_len) % segment_inf[i].sb.size();
                        seg_details[cur_track].sb[j] = segment_inf[i].sb[index];
                    }
                }

                seg_details[cur_track].Rmetal = segment_inf[i].Rmetal;
                seg_details[cur_track].Cmetal = segment_inf[i].Cmetal;


                if (0 == ipart) {
                    seg_details[cur_track].arch_wire_switch = arch_wire_switch;
                    seg_details[cur_track].arch_opin_switch = arch_opin_switch;
                    seg_details[cur_track].first_opin_switch = arch_opin_switch;
                    seg_details[cur_track].first_wire_switch = arch_wire_switch;
                    first_part_opin_switch = arch_opin_switch;
                    first_part_wire_switch = arch_wire_switch;
                } else {
                    seg_details[cur_track].arch_wire_switch = -1;
                    seg_details[cur_track].arch_opin_switch = -1;
                    seg_details[cur_track].first_wire_switch = first_part_wire_switch;
                    seg_details[cur_track].first_opin_switch = first_part_opin_switch;
                }

                if (BI_DIRECTIONAL == directionality) {
                    seg_details[cur_track].direction = BI_DIRECTION;
                } else {
                    VTR_ASSERT(UNI_DIRECTIONAL == directionality);
                    seg_details[cur_track].direction = (itrack % 2) ? DEC_DIRECTION : INC_DIRECTION;
                }

                seg_details[cur_track].index = i;
                seg_details[cur_track].cost_index = i + index_offset + ipart;
                
                ++cur_track;
            }
            tmp_sum_len += tmp_seg_len[ipart];
            cost_index_num = i + index_offset + ipart + 1;
        }
        if (tmp_seg_len.size() > 1)
            index_offset += tmp_seg_len.size() - 1;
    } /* End for each segment type. */

    for(int i = 0; i < 1; i++){
        first_track = cur_track;
        ntracks = medium_num_xory;
        if (ntracks < 1) {
            continue;
        }
        /* Avoid divide by 0 if ntracks */
        longline = false;
        length = 1;

        arch_wire_switch = arch_switch;
        arch_opin_switch = arch_switch;
        VTR_ASSERT((arch_wire_switch == arch_opin_switch) || (directionality != UNI_DIRECTIONAL));

        is_gsb_medium = true;

        medium_seg_name=medium_seg_names;

        /* Set up the tracks of same type */
        group_start = 0;
        for (itrack = 0; itrack < ntracks; itrack++) {
            /* set the name of the segment type this track belongs to */
            seg_details[cur_track].type_name = medium_seg_name;
            seg_details[cur_track].part_idx = 0;
            /* Remember the start track of the current wire group */
            if ((itrack / fac) % length == 0 && (itrack % fac) == 0) {
                group_start = cur_track;
            }

            seg_details[cur_track].length = 1;
            seg_details[cur_track].longline = false;
            seg_details[cur_track].switch_dir_type = NORMAL_TYPE;
            seg_details[cur_track].bend_len = 1;

            seg_details[cur_track].is_gsb_medium = is_gsb_medium;

            /* Stagger the start points in for each track set. The
             * pin mappings should be aware of this when chosing an
             * intelligent way of connecting pins and tracks.
             * cur_track is used as an offset so that extra tracks
             * from different segment types are hopefully better
             * balanced. */
            seg_details[cur_track].start = 1;

            /* These properties are used for vpr_to_phy_track to determine
             * * twisting of wires. */
            seg_details[cur_track].group_start = group_start;
            seg_details[cur_track].group_size = fac;
            seg_details[cur_track].seg_start = -1;
            seg_details[cur_track].seg_end = -1;
            seg_details[cur_track].Rmetal = 0;
            seg_details[cur_track].Cmetal = 0;


            seg_details[cur_track].arch_wire_switch = arch_wire_switch;
            seg_details[cur_track].arch_opin_switch = arch_opin_switch;
            seg_details[cur_track].first_opin_switch = arch_opin_switch;
            seg_details[cur_track].first_wire_switch = arch_wire_switch;

            if (BI_DIRECTIONAL == directionality) {
                seg_details[cur_track].direction = BI_DIRECTION;
            } else {
                VTR_ASSERT(UNI_DIRECTIONAL == directionality);
                seg_details[cur_track].direction = (itrack % 2) ? DEC_DIRECTION : INC_DIRECTION;
            }

            seg_details[cur_track].index = segment_inf.size() -1 + i;
            seg_details[cur_track].cost_index = cost_index_num + i;
            VTR_ASSERT(0 == strcmp(segment_inf[segment_inf.size() -1 + i].name.c_str(), medium_seg_name.c_str()));

            ++cur_track;
        }
    }

    /* free variables */
    vtr::free(sets_per_seg_type);

    if (num_seg_details) {
        *num_seg_details = cur_track;
    }
    return seg_details;
}

/* Allocates and loads the chan_details data structure, a 2D array of
 * seg_details structures. This array is used to handle unique seg_details
 * (ie. channel segments) for each horizontal and vertical channel. */

void alloc_and_load_chan_details(const DeviceGrid& grid,
                                 const t_chan_width* nodes_per_chan,
                                 const bool trim_empty_channels,
                                 const bool trim_obs_channels,
                                 const int num_seg_details,
                                 const t_seg_details* seg_details,
                                 const t_medium_seg_inf medium_seg_inf,
                                 t_chan_details& chan_details_x,
                                 t_chan_details& chan_details_y) {
    chan_details_x = init_chan_details(grid, nodes_per_chan,
                                       num_seg_details, seg_details, medium_seg_inf, SEG_DETAILS_X);
    chan_details_y = init_chan_details(grid, nodes_per_chan,
                                       num_seg_details, seg_details, medium_seg_inf, SEG_DETAILS_Y);

    /* Obstruct channel segment details based on grid block widths/heights */
    obstruct_chan_details(grid, nodes_per_chan, medium_seg_inf,
                          trim_empty_channels, trim_obs_channels,
                          chan_details_x, chan_details_y);

    /* Adjust segment start/end based on obstructed channels, if any */
    adjust_chan_details(grid, nodes_per_chan,
                        chan_details_x, chan_details_y);
}

t_chan_details init_chan_details(const DeviceGrid& grid,
                                 const t_chan_width* nodes_per_chan,
                                 const int num_seg_details,
                                 const t_seg_details* seg_details,
                                 const t_medium_seg_inf medium_seg_inf,
                                 const enum e_seg_details_type seg_details_type) {
    auto& device_ctx = g_vpr_ctx.device();
    int medium_num_xory = medium_seg_inf.gsb_mdeium_num_xory;
    VTR_ASSERT(num_seg_details <= nodes_per_chan->max + medium_num_xory);

    t_chan_details chan_details({grid.width(), grid.height(), size_t(num_seg_details)});

    for (size_t x = 0; x < grid.width(); ++x) {
        for (size_t y = 0; y < grid.height(); ++y) {
            t_chan_seg_details* p_seg_details = chan_details[x][y].data();
            for (int i = 0; i < num_seg_details; ++i) {
                p_seg_details[i] = t_chan_seg_details(&seg_details[i]);
                
                int seg_start = -1;
                int seg_end = -1;

                if (seg_details_type == SEG_DETAILS_X) {
                    seg_start = get_seg_start(p_seg_details, i, y, x);
                    seg_end = get_seg_end(p_seg_details, i, seg_start, y, grid.width() - 2); //-2 for no perim channels
                }
                if (seg_details_type == SEG_DETAILS_Y) {
                    seg_start = get_seg_start(p_seg_details, i, x, y);
                    seg_end = get_seg_end(p_seg_details, i, seg_start, x, grid.height() - 2); //-2 for no perim channels
                }

                p_seg_details[i].set_seg_start(seg_start);
                p_seg_details[i].set_seg_end(seg_end);

                if (seg_details_type == SEG_DETAILS_X) {
                    if(i >= nodes_per_chan->x_list[y] && i < num_seg_details - medium_num_xory) {
                        p_seg_details[i].set_length(0);
                    }
                }
                if (seg_details_type == SEG_DETAILS_Y) {
                    if (i >= nodes_per_chan->y_list[x] && i < num_seg_details - medium_num_xory) {
                        p_seg_details[i].set_length(0);
                    }
                }
            }
        }
    }
    return chan_details;
}

void obstruct_chan_details(const DeviceGrid& grid,
                           const t_chan_width* nodes_per_chan,
                           const t_medium_seg_inf medium_seg_inf,
                           const bool trim_empty_channels,
                           const bool trim_obs_channels,
                           t_chan_details& chan_details_x,
                           t_chan_details& chan_details_y) {
    auto& device_ctx = g_vpr_ctx.device();
    int track_max = nodes_per_chan->max + medium_seg_inf.gsb_mdeium_num_xory;

    /* Iterate grid to find and obstruct based on multi-width/height blocks */
    for (size_t x = 0; x < grid.width() - 1; ++x) {
        for (size_t y = 0; y < grid.height() - 1; ++y) {
            if (!trim_obs_channels)
                continue;

            if (grid[x][y].type == device_ctx.EMPTY_TYPE)
                continue;
            if (grid[x][y].width_offset > 0 || grid[x][y].height_offset > 0)
                continue;
            if (grid[x][y].type->width == 1 && grid[x][y].type->height == 1)
                continue;

            if (grid[x][y].type->height > 1) {
                for (int dx = 0; dx <= grid[x][y].type->width - 1; ++dx) {
                    for (int dy = 0; dy < grid[x][y].type->height - 1; ++dy) {
                        for (int track = 0; track < track_max; ++track) {
                            chan_details_x[x + dx][y + dy][track].set_length(0);
                        }
                    }
                }
            }
            if (grid[x][y].type->width > 1) {
                for (int dy = 0; dy <= grid[x][y].type->height - 1; ++dy) {
                    for (int dx = 0; dx < grid[x][y].type->width - 1; ++dx) {
                        for (int track = 0; track < track_max; ++track) {
                            chan_details_y[x + dx][y + dy][track].set_length(0);
                        }
                    }
                }
            }
        }
    }

    /* Iterate grid again to find and obstruct based on neighboring EMPTY and/or IO types */
    for (size_t x = 0; x <= grid.width() - 2; ++x) {      //-2 for no perim channels
        for (size_t y = 0; y <= grid.height() - 2; ++y) { //-2 for no perim channels

            if (!trim_empty_channels)
                continue;

            if (is_io_type(grid[x][y].type)) {
                if ((x == 0) || (y == 0))
                    continue;
            }
            if (grid[x][y].type == device_ctx.EMPTY_TYPE) {
                if ((x == grid.width() - 2) && is_io_type(grid[x + 1][y].type)) //-2 for no perim channels
                    continue;
                if ((y == grid.height() - 2) && is_io_type(grid[x][y + 1].type)) //-2 for no perim channels
                    continue;
            }

            if (is_io_type(grid[x][y].type) || (grid[x][y].type == device_ctx.EMPTY_TYPE)) {
                if (is_io_type(grid[x][y + 1].type) || (grid[x][y + 1].type == device_ctx.EMPTY_TYPE)) {
                    for (int track = 0; track < track_max; ++track) {
                        chan_details_x[x][y][track].set_length(0);
                    }
                }
            }
            if (is_io_type(grid[x][y].type) || (grid[x][y].type == device_ctx.EMPTY_TYPE)) {
                if (is_io_type(grid[x + 1][y].type) || (grid[x + 1][y].type == device_ctx.EMPTY_TYPE)) {
                    for (int track = 0; track < track_max; ++track) {
                        chan_details_y[x][y][track].set_length(0);
                    }
                }
            }
        }
    }
}

void adjust_chan_details(const DeviceGrid& grid,
                         const t_chan_width* nodes_per_chan,
                         t_chan_details& chan_details_x,
                         t_chan_details& chan_details_y) {
    for (size_t y = 0; y <= grid.height() - 2; ++y) {    //-2 for no perim channels
        for (size_t x = 0; x <= grid.width() - 2; ++x) { //-2 for no perim channels

            /* Ignore any non-obstructed channel seg_detail structures */
            if (chan_details_x[x][y][0].length() > 0)
                continue;

            adjust_seg_details(x, y, grid, nodes_per_chan,
                               chan_details_x, SEG_DETAILS_X);
        }
    }

    for (size_t x = 0; x <= grid.width() - 2; ++x) {      //-2 for no perim channels
        for (size_t y = 0; y <= grid.height() - 2; ++y) { //-2 for no perim channels

            /* Ignore any non-obstructed channel seg_detail structures */
            if (chan_details_y[x][y][0].length() > 0)
                continue;

            adjust_seg_details(x, y, grid, nodes_per_chan,
                               chan_details_y, SEG_DETAILS_Y);
        }
    }
}

void adjust_seg_details(const int x,
                        const int y,
                        const DeviceGrid& grid,
                        const t_chan_width* nodes_per_chan,
                        t_chan_details& chan_details,
                        const enum e_seg_details_type seg_details_type) {
    int seg_index = (seg_details_type == SEG_DETAILS_X ? x : y);

    for (int track = 0; track < nodes_per_chan->max; ++track) {
        int lx = (seg_details_type == SEG_DETAILS_X ? x - 1 : x);
        int ly = (seg_details_type == SEG_DETAILS_X ? y : y - 1);
        if (lx < 0 || ly < 0 || chan_details[lx][ly][track].length() == 0)
            continue;

        while (chan_details[lx][ly][track].seg_end() >= seg_index) {
            chan_details[lx][ly][track].set_seg_end(seg_index - 1);
            lx = (seg_details_type == SEG_DETAILS_X ? lx - 1 : lx);
            ly = (seg_details_type == SEG_DETAILS_X ? ly : ly - 1);
            if (lx < 0 || ly < 0 || chan_details[lx][ly][track].length() == 0)
                break;
        }
    }

    for (int track = 0; track < nodes_per_chan->max; ++track) {
        size_t lx = (seg_details_type == SEG_DETAILS_X ? x + 1 : x);
        size_t ly = (seg_details_type == SEG_DETAILS_X ? y : y + 1);
        if (lx > grid.width() - 2 || ly > grid.height() - 2 || chan_details[lx][ly][track].length() == 0) //-2 for no perim channels
            continue;

        while (chan_details[lx][ly][track].seg_start() <= seg_index) {
            chan_details[lx][ly][track].set_seg_start(seg_index + 1);
            lx = (seg_details_type == SEG_DETAILS_X ? lx + 1 : lx);
            ly = (seg_details_type == SEG_DETAILS_X ? ly : ly + 1);
            if (lx > grid.width() - 2 || ly > grid.height() - 2 || chan_details[lx][ly][track].length() == 0) //-2 for no perim channels
                break;
        }
    }
}

void free_chan_details(t_chan_details& chan_details_x,
                       t_chan_details& chan_details_y) {
    chan_details_x.clear();
    chan_details_y.clear();
}

/* Returns the segment number at which the segment this track lies on        *
 * started.                                                                  */
int get_seg_start(const t_chan_seg_details* seg_details,
                  const int itrack,
                  const int chan_num,
                  const int seg_num) {
    int seg_start = 0;
    if (seg_details[itrack].seg_start() >= 0) {
        seg_start = seg_details[itrack].seg_start();

    } else {
        seg_start = 1;
        if (false == seg_details[itrack].longline()) {
            int length = seg_details[itrack].length();
            int start = seg_details[itrack].start();

            /* Start is guaranteed to be between 1 and length.  Hence adding length to *
             * the quantity in brackets below guarantees it will be nonnegative.       */

            VTR_ASSERT(start > 0);
            VTR_ASSERT(start <= length);

            /* NOTE: Start points are staggered between different channels.
             * The start point must stagger backwards as chan_num increases.
             * Unidirectional routing expects this to allow the N-to-N
             * assumption to be made with respect to ending wires in the core. */
            seg_start = seg_num - (seg_num + length + chan_num - start) % length;
            if (seg_start < 1) {
                seg_start = 1;
            }
        }
    }
    return seg_start;
}

int get_seg_end(const t_chan_seg_details* seg_details, const int itrack, const int istart, const int chan_num, const int seg_max) {
    if (seg_details[itrack].longline()) {
        return seg_max;
    }

    if (seg_details[itrack].seg_end() >= 0) {
        return seg_details[itrack].seg_end();
    }

    int len = seg_details[itrack].length();
    int ofs = seg_details[itrack].start();

    /* Normal endpoint */
    int seg_end = istart + len - 1;

    /* If start is against edge it may have been clipped */
    if (1 == istart) {
        /* If the (staggered) startpoint of first full wire wasn't
         * also 1, we must be the clipped wire */
        int first_full = (len - (chan_num % len) + ofs - 1) % len + 1;
        if (first_full > 1) {
            /* then we stop just before the first full seg */
            seg_end = first_full - 1;
        }
    }

    /* Clip against far edge */
    if (seg_end > seg_max) {
        seg_end = seg_max;
    }

    return seg_end;
}

/* Returns the number of tracks to which clb opin #ipin at (i,j) connects.   *
 * Also stores the nodes to which this pin connects in rr_edges_to_create    */
int get_bidir_opin_connections(const int i,
                               const int j,
                               const int ipin,
                               const int from_rr_node,
                               t_rr_edge_info_set& rr_edges_to_create,
                               const t_pin_to_track_lookup_bi& opin_to_track_map,
                               const t_rr_node_indices& L_rr_node_indices,
                               const t_chan_details& chan_details_x,
                               const t_chan_details& chan_details_y) {
    int num_conn, tr_i, tr_j, chan, seg;
    int to_switch, to_node;
    int is_connected_track;
    t_type_ptr type;
    t_rr_type to_type;

    auto& device_ctx = g_vpr_ctx.device();

    type = device_ctx.grid[i][j].type;
    int width_offset = device_ctx.grid[i][j].width_offset;
    int height_offset = device_ctx.grid[i][j].height_offset;

    num_conn = 0;

    /* [0..device_ctx.num_block_types-1][0..num_pins-1][0..width][0..height][0..3][0..Fc-1] */
    for (e_side side : SIDES) {
        /* Figure out coords of channel segment based on side */
        tr_i = ((side == LEFT) ? (i - 1) : i);
        tr_j = ((side == BOTTOM) ? (j - 1) : j);

        to_type = ((side == LEFT) || (side == RIGHT)) ? CHANY : CHANX;

        chan = ((to_type == CHANX) ? tr_j : tr_i);
        seg = ((to_type == CHANX) ? tr_i : tr_j);

        bool vert = !((side == TOP) || (side == BOTTOM));

        /* Don't connect where no tracks on fringes */
        if ((tr_i < 0) || (tr_i > int(device_ctx.grid.width() - 2))) { //-2 for no perimeter channels
            continue;
        }
        if ((tr_j < 0) || (tr_j > int(device_ctx.grid.height() - 2))) { //-2 for no perimeter channels
            continue;
        }
        if ((CHANX == to_type) && (tr_i < 1)) {
            continue;
        }
        if ((CHANY == to_type) && (tr_j < 1)) {
            continue;
        }
        if (opin_to_track_map[type->index].empty()) {
            continue;
        }

        is_connected_track = false;

        const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg] : chan_details_x[seg][chan]).data();

        /* Iterate of the opin to track connections */
        for (int to_track : opin_to_track_map[type->index][ipin][width_offset][height_offset][side]) {
            /* Skip unconnected connections */
            if (OPEN == to_track || is_connected_track) {
                is_connected_track = true;
                VTR_ASSERT(OPEN == opin_to_track_map[type->index][ipin][width_offset][height_offset][side][0]);
                continue;
            }

            /* Only connect to wire if there is a CB */
            if (is_cblock(chan, seg, to_track, seg_details)) {
                to_switch = seg_details[to_track].arch_wire_switch();
                to_node = get_rr_node_index(L_rr_node_indices, tr_i, tr_j, to_type, to_track);

                if (to_node == OPEN) {
                    continue;
                }

                rr_edges_to_create.emplace_back(from_rr_node, to_node, to_switch);
                ++num_conn;
            }
        }
    }

    return num_conn;
}

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
                                    const DeviceGrid& grid) {
    /* Gets a linked list of Fc nodes of specified seg_type_index to connect
     * to in given chan seg. Fc_ofs is used for the opin staggering pattern. */
    //modified by shikc in 200220
    int* inc_muxes = nullptr;
    //int* dec_muxes = nullptr;
    int num_inc_muxes,  iconn;
    int inc_inode_index;
    int inc_mux;
    int inc_track;
    int x, y;
    int num_edges;

    *Fc_clipped = false;

    /* Fc is assigned in pairs so check it is even. */
    VTR_ASSERT(Fc % 2 == 0);

    /* get_rr_node_indices needs x and y coords. */
    x = ((CHANX == chan_type) ? seg : chan);
    y = ((CHANX == chan_type) ? chan : seg);
    //modified by shikc in 200221
    /* Get the lists of possible muxes. */
    int dummy;
    e_side inc_side = (CHANX == chan_type ? RIGHT : TOP);
    //e_side dec_side = (CHANX == chan_type ? LEFT : BOTTOM);
    bool is_dangling_side_inc = label_dangling_side(inc_side, x, y, grid);
    //bool is_dangling_side_dec = label_dangling_side(dec_side, x, y, grid);
    inc_muxes = label_wire_muxes(chan, seg, seg_details, seg_type_index, max_len,
    							  INC_DIRECTION, max_chan_width, true, &num_inc_muxes, &dummy, is_dangling_side_inc);
    //dec_muxes = label_wire_muxes(chan, seg - 1, seg_details, seg_type_index, max_len,
    //                             DEC_DIRECTION, max_chan_width, true, &num_dec_muxes, &dummy)

    /* Clip Fc to the number of muxes. */
    if (((Fc / 2) > num_inc_muxes)) {
        *Fc_clipped = true;
        Fc = 2 * num_inc_muxes;
    }

    /* Assign tracks to meet Fc demand */
    num_edges = 0;
    //Fc_ofs[chan][seg][seg_type_index][0] = 0;
    for (iconn = 0; iconn < (Fc / 2); ++iconn) {
        /* Figure of the next mux to use for the 'inc' and 'dec' connections */
        inc_mux = Fc_ofs[chan][seg][seg_type_index][0] % num_inc_muxes;
        //dec_mux = Fc_ofs[chan][seg][seg_type_index] % num_dec_muxes;
        ++Fc_ofs[chan][seg][seg_type_index][0];

        /* Figure out the track it corresponds to. */
        VTR_ASSERT(inc_muxes != nullptr);
        inc_track = inc_muxes[inc_mux];

        //VTR_ASSERT(dec_muxes != nullptr);
        //dec_track = dec_muxes[dec_mux];

        /* Figure the inodes of those muxes */
        inc_inode_index = get_rr_node_index(L_rr_node_indices, x, y, chan_type, inc_track);
        //dec_inode_index = get_rr_node_index(L_rr_node_indices, x - 1, y, chan_type, dec_track);

        if (inc_inode_index == OPEN) {
            continue;
        }

        /* Add to the list. */
        if (seg_details[inc_track].part_idx() > 0 && is_dangling_side_inc)
            rr_edges_to_create.emplace_back(from_rr_node, inc_inode_index, seg_details[inc_track].first_opin_switch());
        else
            rr_edges_to_create.emplace_back(from_rr_node, inc_inode_index, seg_details[inc_track].arch_opin_switch());
        ++num_edges;

        //rr_edges_to_create.emplace_back(from_rr_node, dec_inode_index, seg_details[dec_track].arch_opin_switch());
        //++num_edges;
    }

    if (inc_muxes) {
        vtr::free(inc_muxes);
        inc_muxes = nullptr;
    }

    return num_edges;
}
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
                                    const DeviceGrid& grid) {
    /* Gets a linked list of Fc nodes of specified seg_type_index to connect
     * to in given chan seg. Fc_ofs is used for the opin staggering pattern. */
    //modified by shikc in 200220
    //int* inc_muxes = nullptr;
    int* dec_muxes = nullptr;
    int num_dec_muxes, iconn;
    int dec_inode_index;
    int  dec_mux;
    int dec_track;
    int x, y;
    int num_edges;

    *Fc_clipped = false;

    /* Fc is assigned in pairs so check it is even. */
    VTR_ASSERT(Fc % 2 == 0);

    /* get_rr_node_indices needs x and y coords. */
    x = ((CHANX == chan_type) ? seg : chan);
    y = ((CHANX == chan_type) ? chan : seg);
    //modified by shikc in 200221
    /* Get the lists of possible muxes. */
    int dummy;
    //e_side inc_side = (CHANX == chan_type ? RIGHT : TOP);
    e_side dec_side = (CHANX == chan_type ? LEFT : BOTTOM);
    //bool is_dangling_side_inc = label_dangling_side(inc_side, x, y, grid);
    bool is_dangling_side_dec = label_dangling_side(dec_side, x, y, grid);
    dec_muxes = label_wire_muxes(chan, seg - 1, seg_details, seg_type_index, max_len,
                                 DEC_DIRECTION, max_chan_width, true, &num_dec_muxes, &dummy, is_dangling_side_dec);

    /* Clip Fc to the number of muxes. */
    if (((Fc / 2) > num_dec_muxes)) {
        *Fc_clipped = true;
        Fc = 2 * num_dec_muxes;
    }

    /* Assign tracks to meet Fc demand */
    num_edges = 0;
    //Fc_ofs[chan][seg][seg_type_index] = 0;
    for (iconn = 0; iconn < (Fc / 2); ++iconn) {
        /* Figure of the next mux to use for the 'inc' and 'dec' connections */
        //inc_mux = Fc_ofs[chan][seg][seg_type_index] % num_inc_muxes;
        dec_mux = Fc_ofs[chan][seg-1][seg_type_index][1] % num_dec_muxes;
        ++Fc_ofs[chan][seg-1][seg_type_index][1];

        /* Figure out the track it corresponds to. */
        //VTR_ASSERT(inc_muxes != nullptr);
        //inc_track = inc_muxes[inc_mux];

        VTR_ASSERT(dec_muxes != nullptr);
        dec_track = dec_muxes[dec_mux];

        /* Figure the inodes of those muxes */
        if(CHANX == chan_type) {
            //inc_inode_index = get_rr_node_index(L_rr_node_indices, x, y, chan_type, inc_track);
            dec_inode_index = get_rr_node_index(L_rr_node_indices, x - 1, y, chan_type, dec_track);
        }
        else{
            //inc_inode_index = get_rr_node_index(L_rr_node_indices, x, y, chan_type, inc_track);
            dec_inode_index = get_rr_node_index(L_rr_node_indices, x, y-1, chan_type, dec_track);}

        if (dec_inode_index == OPEN) {
            continue;
        }

        /* Add to the list. */
        //rr_edges_to_create.emplace_back(from_rr_node, inc_inode_index, seg_details[inc_track].arch_opin_switch());
        //++num_edges;
        if (seg_details[dec_track].part_idx() > 0 && is_dangling_side_dec)
            rr_edges_to_create.emplace_back(from_rr_node, dec_inode_index, seg_details[dec_track].first_opin_switch());
        else
            rr_edges_to_create.emplace_back(from_rr_node, dec_inode_index, seg_details[dec_track].arch_opin_switch());
        ++num_edges;
    }


    if (dec_muxes) {
        vtr::free(dec_muxes);
        dec_muxes = nullptr;
    }

    return num_edges;
}

bool is_cblock(const int chan, const int seg, const int track, const t_chan_seg_details* seg_details) {
    int length, ofs, start_seg;

    length = seg_details[track].length();

    /* Make sure they gave us correct start */
    start_seg = get_seg_start(seg_details, track, chan, seg);

    ofs = seg - start_seg;

    VTR_ASSERT(ofs >= 0);
    VTR_ASSERT(ofs < length);

    /* If unidir segment that is going backwards, we need to flip the ofs */
    if (DEC_DIRECTION == seg_details[track].direction()) {
        ofs = (length - 1) - ofs;
    }

    return seg_details[track].cb(ofs);
}

void dump_seg_details(const t_chan_seg_details* seg_details,
                      int max_chan_width,
                      FILE* fp) {
    for (int i = 0; i < max_chan_width; i++) {
        fprintf(fp, "track: %d\n", i);
        fprintf(fp, "length: %d  start: %d",
                seg_details[i].length(), seg_details[i].start());

        if (seg_details[i].length() > 0) {
            if (seg_details[i].seg_start() >= 0 && seg_details[i].seg_end() >= 0) {
                fprintf(fp, " [%d,%d]",
                        seg_details[i].seg_start(), seg_details[i].seg_end());
            }
            fprintf(fp, "  longline: %d  arch_wire_switch: %d  arch_opin_switch: %d",
                    seg_details[i].longline(),
                    seg_details[i].arch_wire_switch(), seg_details[i].arch_opin_switch());
        }
        fprintf(fp, "\n");

        fprintf(fp, "Rmetal: %g  Cmetal: %g\n",
                seg_details[i].Rmetal(), seg_details[i].Cmetal());

        fprintf(fp, "direction: %s\n",
                DIRECTION_STRING[seg_details[i].direction()]);

        fprintf(fp, "name: %s\n", seg_details[i].type_name().c_str());

        fprintf(fp, "cb list:  ");
        for (int j = 0; j < seg_details[i].length(); j++)
            fprintf(fp, "%d ", seg_details[i].cb(j));
        fprintf(fp, "\n");

        fprintf(fp, "sb list: ");
        for (int j = 0; j <= seg_details[i].length(); j++)
            fprintf(fp, "%d ", seg_details[i].sb(j));
        fprintf(fp, "\n");

        fprintf(fp, "bend part index: %d (max: %d)\n", seg_details[i].part_idx(), seg_details[i].bend_len());
        fprintf(fp, "seg index: %d\n", seg_details[i].index());
        fprintf(fp, "cost index: %d\n", seg_details[i].cost_index());

        if (seg_details[i].switch_dir_type() == UP_TYPE)
            fprintf(fp, "bend switch type : UP ");
        else if (seg_details[i].switch_dir_type() == DOWN_TYPE)
            fprintf(fp, "bend switch type : DOWN ");
        else if (seg_details[i].switch_dir_type() == NORMAL_TYPE)
            fprintf(fp, "bend switch type : NORMAL_TYPE ");

        fprintf(fp, "first_wire_switch: %d first_opin_switch: %d\n", seg_details[i].first_wire_switch(), seg_details[i].first_opin_switch());

        fprintf(fp, "\n");
    }
}

/* Dumps out an array of seg_details structures to file fname.  Used only   *
 * for debugging.                                                           */
void dump_seg_details(const t_chan_seg_details* seg_details,
                      int max_chan_width,
                      const char* fname) {
    FILE* fp = vtr::fopen(fname, "w");
    dump_seg_details(seg_details, max_chan_width, fp);
    fclose(fp);
}

/* Dumps out a 2D array of chan_details structures to file fname.  Used     *
 * only for debugging.                                                      */
void dump_chan_details(const t_chan_details& chan_details_x,
                       const t_chan_details& chan_details_y,
                       int max_chan_width,
                       const DeviceGrid& grid,
                       const char* fname) {
    FILE* fp = vtr::fopen(fname, "w");
    if (fp) {
        for (size_t y = 0; y <= grid.height() - 2; ++y) {    //-2 for no perim channels
            for (size_t x = 0; x <= grid.width() - 2; ++x) { //-2 for no perim channels

                fprintf(fp, "========================\n");
                fprintf(fp, "chan_details_x: [%zu][%zu]\n", x, y);
                fprintf(fp, "========================\n");

                const t_chan_seg_details* seg_details = chan_details_x[x][y].data();
                dump_seg_details(seg_details, max_chan_width, fp);
            }
        }
        for (size_t x = 0; x <= grid.width() - 2; ++x) {      //-2 for no perim channels
            for (size_t y = 0; y <= grid.height() - 2; ++y) { //-2 for no perim channels

                fprintf(fp, "========================\n");
                fprintf(fp, "chan_details_y: [%zu][%zu]\n", x, y);
                fprintf(fp, "========================\n");

                const t_chan_seg_details* seg_details = chan_details_y[x][y].data();
                dump_seg_details(seg_details, max_chan_width, fp);
            }
        }
    }
    fclose(fp);
}

/* Dumps out a 2D array of switch block pattern structures to file fname. *
 * Used for debugging purposes only.                                      */
void dump_sblock_pattern(const t_sblock_pattern& sblock_pattern,
                         int max_chan_width,
                         const DeviceGrid& grid,
                         const char* fname) {
    FILE* fp = vtr::fopen(fname, "w");
    if (fp) {
        for (size_t y = 0; y <= grid.height() - 2; ++y) {
            for (size_t x = 0; x <= grid.width() - 2; ++x) {
                fprintf(fp, "==========================\n");
                fprintf(fp, "sblock_pattern: [%zu][%zu]\n", x, y);
                fprintf(fp, "==========================\n");

                for (int from_side = 0; from_side < 4; ++from_side) {
                    for (int to_side = 0; to_side < 4; ++to_side) {
                        if (from_side == to_side)
                            continue;

                        const char* psz_from_side = "?";
                        switch (from_side) {
                            case 0:
                                psz_from_side = "T";
                                break;
                            case 1:
                                psz_from_side = "R";
                                break;
                            case 2:
                                psz_from_side = "B";
                                break;
                            case 3:
                                psz_from_side = "L";
                                break;
                            default:
                                VTR_ASSERT_MSG(false, "Unrecognized from side");
                                break;
                        }
                        const char* psz_to_side = "?";
                        switch (to_side) {
                            case 0:
                                psz_to_side = "T";
                                break;
                            case 1:
                                psz_to_side = "R";
                                break;
                            case 2:
                                psz_to_side = "B";
                                break;
                            case 3:
                                psz_to_side = "L";
                                break;
                            default:
                                VTR_ASSERT_MSG(false, "Unrecognized to side");
                                break;
                        }

                        for (int from_track = 0; from_track < max_chan_width; ++from_track) {
                            short to_mux = sblock_pattern[x][y][from_side][to_side][from_track][0];
                            short to_track = sblock_pattern[x][y][from_side][to_side][from_track][1];
                            short alt_mux = sblock_pattern[x][y][from_side][to_side][from_track][2];
                            short alt_track = sblock_pattern[x][y][from_side][to_side][from_track][3];

                            if (to_mux == UN_SET && to_track == UN_SET)
                                continue;

                            if (alt_mux == UN_SET && alt_track == UN_SET) {
                                fprintf(fp, "%s %d => %s [%d][%d]\n",
                                        psz_from_side, from_track, psz_to_side,
                                        to_mux, to_track);
                            } else {
                                fprintf(fp, "%s %d => %s [%d][%d] [%d][%d]\n",
                                        psz_from_side, from_track, psz_to_side,
                                        to_mux, to_track, alt_mux, alt_track);
                            }
                        }
                    }
                }
            }
        }
    }
    fclose(fp);
}

void dump_custom_sblock_pattern(
    t_sb_connection_map* sb_conn_map,
    int max_chan_width,
    const DeviceGrid& grid,
    const char* fname) {
    FILE* fp = vtr::fopen(fname, "w");
    if (fp) {
        for (size_t y = 0; y <= grid.height() - 2; ++y) {
            for (size_t x = 0; x <= grid.width() - 2; ++x) {
                fprintf(fp, "==========================\n");
                fprintf(fp, "sblock_pattern: [%zu][%zu]\n", x, y);
                fprintf(fp, "==========================\n");

                for (e_side from_side : SIDES) {
                    for (e_side to_side : SIDES) {
                        if (from_side == to_side)
                            continue;

                        const char* psz_from_side = "?";
                        switch (from_side) {
                            case TOP:
                                psz_from_side = "T";
                                break;
                            case RIGHT:
                                psz_from_side = "R";
                                break;
                            case BOTTOM:
                                psz_from_side = "B";
                                break;
                            case LEFT:
                                psz_from_side = "L";
                                break;
                            default:
                                VTR_ASSERT_MSG(false, "Unrecognized from side");
                                break;
                        }
                        const char* psz_to_side = "?";
                        switch (to_side) {
                            case TOP:
                                psz_to_side = "T";
                                break;
                            case RIGHT:
                                psz_to_side = "R";
                                break;
                            case BOTTOM:
                                psz_to_side = "B";
                                break;
                            case LEFT:
                                psz_to_side = "L";
                                break;
                            default:
                                VTR_ASSERT_MSG(false, "Unrecognized to side");
                                break;
                        }

                        Switchblock_Lookup sb_conn(x, y, from_side, to_side);
                        for (auto edge : (*sb_conn_map)[sb_conn]) {
                            fprintf(fp, "%s %hd => %s %hd  [switch: %hd]\n",
                                    psz_from_side, edge.from_wire, psz_to_side,
                                    edge.to_wire, edge.switch_ind);
                        }
                    }
                }
            }
        }
    }
    fclose(fp);
}

static void load_chan_rr_indices(const int max_chan_width,
                                 const int chan_len,
                                 const int num_chans,
                                 const t_rr_type type,
                                 const t_chan_details& chan_details,
                                 t_rr_node_indices& indices,
                                 const DeviceGrid& grid,
                                 int* index) {
    VTR_ASSERT(indices[type].size() == size_t(num_chans));
    auto& device_ctx = g_vpr_ctx.device();
    for (int chan = 0; chan < num_chans - 1; ++chan) {
        VTR_ASSERT(indices[type][chan].size() == size_t(chan_len));

        for (int seg = 1; seg < chan_len - 1; ++seg) {
            VTR_ASSERT(indices[type][chan][seg].size() == NUM_SIDES);

            
            indices[type][chan][seg][0].resize(max_chan_width, OPEN);
        }
    }

    for (int chan = 0; chan < num_chans - 1; ++chan) {
        for (int seg = 1; seg < chan_len - 1; ++seg) {
            /* Assign an inode to the starts of tracks */
            int x = (type == CHANX ? seg : chan);
            int y = (type == CHANX ? chan : seg);
            if((x == 0 || y == 0) && grid[x][y].type == device_ctx.EMPTY_TYPE){
                continue;
            }
            const t_chan_seg_details* seg_details = chan_details[x][y].data();

            for (unsigned track = 0; track < indices[type][chan][seg][0].size(); ++track) {
                if (seg_details[track].length() <= 0)
                    continue;
                
                if(seg_details[track].is_gsb_medium())
                    continue;

                int start = get_seg_start(seg_details, track, chan, seg);

                /* If the start of the wire doesn't have a inode,
                 * assign one to it. */
                int inode = indices[type][chan][start][0][track];
                if (OPEN == inode) {
                    inode = *index;
                    ++(*index);

                    indices[type][chan][start][0][track] = inode;
                }

                /* Assign inode of start of wire to current position */
                indices[type][chan][seg][0][track] = inode;
            }
        }
    }
}

static void load_block_rr_indices(const DeviceGrid& grid,
                                  t_rr_node_indices& indices,
                                  int* index) {
    //Walk through the grid assigning indices to SOURCE/SINK IPIN/OPIN
    for (size_t x = 0; x < grid.width(); x++) {
        for (size_t y = 0; y < grid.height(); y++) {
            if (grid[x][y].width_offset == 0 && grid[x][y].height_offset == 0) {
                //Process each block from it's root location
                t_type_ptr type = grid[x][y].type;

                //Assign indicies for SINKs and SOURCEs
                // Note that SINKS/SOURCES have no side, so we always use side 0
                for (int iclass = 0; iclass < type->num_class; ++iclass) {
                    auto class_type = type->class_inf[iclass].type;
                    if (class_type == DRIVER) {
                        indices[SOURCE][x][y][0].push_back(*index);
                        indices[SINK][x][y][0].push_back(OPEN);
                    } else {
                        VTR_ASSERT(class_type == RECEIVER);
                        indices[SINK][x][y][0].push_back(*index);
                        indices[SOURCE][x][y][0].push_back(OPEN);
                    }
                    ++(*index);
                }
                VTR_ASSERT(indices[SOURCE][x][y][0].size() == size_t(type->num_class));
                VTR_ASSERT(indices[SINK][x][y][0].size() == size_t(type->num_class));

                //Assign indicies for IPINs and OPINs at all offsets from root
                for (int ipin = 0; ipin < type->num_pins; ++ipin) {
                    for (int width_offset = 0; width_offset < type->width; ++width_offset) {
                        int x_tile = x + width_offset;
                        for (int height_offset = 0; height_offset < type->height; ++height_offset) {
                            int y_tile = y + height_offset;
                            for (e_side side : SIDES) {
                                if (type->pinloc[width_offset][height_offset][side][ipin]) {
                                    int iclass = type->pin_class[ipin];
                                    auto class_type = type->class_inf[iclass].type;

                                    if (class_type == DRIVER) {
                                        indices[OPIN][x_tile][y_tile][side].push_back(*index);
                                        indices[IPIN][x_tile][y_tile][side].push_back(OPEN);
                                    } else {
                                        VTR_ASSERT(class_type == RECEIVER);
                                        indices[IPIN][x_tile][y_tile][side].push_back(*index);
                                        indices[OPIN][x_tile][y_tile][side].push_back(OPEN);
                                    }
                                    ++(*index);
                                } else {
                                    indices[IPIN][x_tile][y_tile][side].push_back(OPEN);
                                    indices[OPIN][x_tile][y_tile][side].push_back(OPEN);
                                }
                            }
                        }
                    }
                }

                //Sanity check
                for (int width_offset = 0; width_offset < type->width; ++width_offset) {
                    int x_tile = x + width_offset;
                    for (int height_offset = 0; height_offset < type->height; ++height_offset) {
                        int y_tile = y + height_offset;
                        for (e_side side : SIDES) {
                            VTR_ASSERT(indices[IPIN][x_tile][y_tile][side].size() == size_t(type->num_pins));
                            VTR_ASSERT(indices[OPIN][x_tile][y_tile][side].size() == size_t(type->num_pins));
                        }
                    }
                }
            }
        }
    }

    //Copy the SOURCE/SINK nodes to all offset positions for blocks with width > 1 and/or height > 1
    // This ensures that look-ups on non-root locations will still find the correct SOURCE/SINK
    for (size_t x = 0; x < grid.width(); x++) {
        for (size_t y = 0; y < grid.height(); y++) {
            int width_offset = grid[x][y].width_offset;
            int height_offset = grid[x][y].height_offset;
            if (width_offset != 0 || height_offset != 0) {
                int root_x = x - width_offset;
                int root_y = y - height_offset;

                indices[SOURCE][x][y] = indices[SOURCE][root_x][root_y];
                indices[SINK][x][y] = indices[SINK][root_x][root_y];
            }
        }
    }
}

t_rr_node_indices alloc_and_load_rr_node_indices(const int max_chan_width,
                                                 const DeviceGrid& grid,
                                                 int* index,
                                                 const t_chan_details& chan_details_x,
                                                 const t_chan_details& chan_details_y) {
    /* Allocates and loads all the structures needed for fast lookups of the   *
     * index of an rr_node.  rr_node_indices is a matrix containing the index  *
     * of the *first* rr_node at a given (i,j) location.                       */

    t_rr_node_indices indices;

    /* Alloc the lookup table */
    indices.resize(NUM_RR_TYPES);
    for (t_rr_type rr_type : RR_TYPES) {
        if (rr_type == CHANX) {
            indices[rr_type].resize(grid.height());
            for (size_t y = 0; y < grid.height(); ++y) {
                indices[rr_type][y].resize(grid.width());
                for (size_t x = 0; x < grid.width(); ++x) {
                    indices[rr_type][y][x].resize(NUM_SIDES);
                }
            }
        } else {
            indices[rr_type].resize(grid.width());
            for (size_t x = 0; x < grid.width(); ++x) {
                indices[rr_type][x].resize(grid.height());
                for (size_t y = 0; y < grid.height(); ++y) {
                    indices[rr_type][x][y].resize(NUM_SIDES);
                }
            }
        }
    }

    /* Assign indices for block nodes */
    load_block_rr_indices(grid, indices, index);

    /* Load the data for x and y channels */
    load_chan_rr_indices(max_chan_width, grid.width(), grid.height(),
                         CHANX, chan_details_x, indices, grid, index);
    load_chan_rr_indices(max_chan_width, grid.height(), grid.width(),
                         CHANY, chan_details_y, indices, grid, index);
    return indices;
}

std::vector<int> get_rr_node_chan_wires_at_location(const t_rr_node_indices& L_rr_node_indices,
                                                    t_rr_type rr_type,
                                                    int x,
                                                    int y) {
    VTR_ASSERT(rr_type == CHANX || rr_type == CHANY);

    /* Currently need to swap x and y for CHANX because of chan, seg convention */
    if (CHANX == rr_type) {
        std::swap(x, y);
    }

    return L_rr_node_indices[rr_type][x][y][SIDES[0]];
}

std::vector<int> get_rr_node_indices(const t_rr_node_indices& L_rr_node_indices,
                                     int x,
                                     int y,
                                     t_rr_type rr_type,
                                     int ptc) {
    /*
     * Like get_rr_node_index() but returns all matching nodes,
     * rather than just the first. This is particularly useful for getting all instances
     * of a specific IPIN/OPIN at a specific gird tile (x,y) location.
     */
    std::vector<int> indices;

    if (rr_type == IPIN || rr_type == OPIN) {
        //For pins we need to look at all the sides of the current grid tile

        for (e_side side : SIDES) {
            int rr_node_index = get_rr_node_index(L_rr_node_indices, x, y, rr_type, ptc, side);

            if (rr_node_index >= 0) {
                indices.push_back(rr_node_index);
            }
        }
    } else {
        //Sides do not effect non-pins so there should only be one per ptc
        int rr_node_index = get_rr_node_index(L_rr_node_indices, x, y, rr_type, ptc);

        if (rr_node_index != OPEN) {
            indices.push_back(rr_node_index);
        }
    }

    return indices;
}

int find_average_rr_node_index(int device_width,
                               int device_height,
                               t_rr_type rr_type,
                               int ptc,
                               const t_rr_node_indices& L_rr_node_indices) {
    /* Find and return the index to a rr_node that is located at the "center" *
     * of the current grid array, if possible.  In the event the "center" of  *
     * the grid array is an EMPTY or IO node, then retry alterate locations.  *
     * Worst case, this function will simply return the 1st non-EMPTY and     *
     * non-IO node.                                                           */

    int inode = get_rr_node_index(L_rr_node_indices, (device_width) / 2, (device_height) / 2,
                                  rr_type, ptc);

    if (inode == OPEN) {
        inode = get_rr_node_index(L_rr_node_indices, (device_width) / 4, (device_height) / 4,
                                  rr_type, ptc);
    }
    if (inode == OPEN) {
        inode = get_rr_node_index(L_rr_node_indices, (device_width) / 4 * 3, (device_height) / 4 * 3,
                                  rr_type, ptc);
    }
    if (inode == OPEN) {
        auto& device_ctx = g_vpr_ctx.device();

        for (int x = 0; x < device_width; ++x) {
            for (int y = 0; y < device_height; ++y) {
                if (device_ctx.grid[x][y].type == device_ctx.EMPTY_TYPE)
                    continue;
                if (is_io_type(device_ctx.grid[x][y].type))
                    continue;

                inode = get_rr_node_index(L_rr_node_indices, x, y, rr_type, ptc);
                if (inode != OPEN)
                    break;
            }
            if (inode != OPEN)
                break;
        }
    }
    return (inode);
}

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
                      enum e_directionality directionality) {
    /*
     * Adds the fan-out edges from wire segment at (chan, seg, track) to adjacent
     * blocks along the wire's length
     */

    int j, pass, iconn, phy_track, end, max_conn, ipin, x, y, num_conn;
    t_type_ptr type;

    auto& device_ctx = g_vpr_ctx.device();

    /* End of this wire */
    end = get_seg_end(seg_details, track, seg, chan, chan_length);

    num_conn = 0;
    //modified by shikc in 200220
    for (j = seg; j <= end; j++) {
        if (is_cblock(chan, j, track, seg_details)) {
            for (pass = 0; pass < 2; ++pass) { //pass == 0 => TOP/RIGHT, pass == 1 => BOTTOM/LEFT
                for (int idx = 0; idx < 4; ++idx) {
                    e_side side;
                    /*if (CHANX == chan_type) {
                    x = j;
                    y = chan + pass;
                    side = (0 == pass ? TOP : BOTTOM);
                } else {
                    VTR_ASSERT(CHANY == chan_type);
                    x = chan + pass;
                    y = j;
                    side = (0 == pass ? RIGHT : LEFT);
                }*/
                    if (CHANX == chan_type) {
                        if(pass == 0) {
                            if (idx == 0) {
                                x = j;
                                y = chan;
                                side = TOP;
                            } else if (idx == 1) {
                                x = j - 1;
                                y = chan + 1;
                                side = BOTTOM;
                            } else if (idx == 2) {
                                x = j - 1;
                                y = chan;
                                side = RIGHT;
                            } else if (idx == 3) {
                                x = j;
                                y = chan + 1;
                                side = LEFT;
                            }
                        } else {
                            if (idx == 0) {
                                x = j;
                                y = chan + 1;
                                side = BOTTOM;
                            } else if (idx == 1) {
                                x = j + 1;
                                y = chan;
                                side = TOP;
                            } else if (idx == 2) {
                                x = j + 1;
                                y = chan + 1;
                                side = LEFT;
                            } else if (idx == 3) {
                                x = j;
                                y = chan;
                                side = RIGHT;
                            }
                        }
                        if(x < 0) continue;

                    } else {
                        VTR_ASSERT(CHANY == chan_type);
                        if (pass == 0) {
                            if (idx == 0) {
                                x = chan;
                                y = j;
                                side = RIGHT;
                            } else if (idx == 1) {
                                x = chan;
                                y = j + 1;
                                side = BOTTOM;
                            } else if (idx == 2) {
                                x = chan + 1;
                                y = j;
                                side = TOP;
                            } else if (idx == 3) {
                                x = chan + 1;
                                y = j + 1;
                                side = LEFT;
                            }
                        } else {
                            if (idx == 0) {
                                x = chan + 1;
                                y = j;
                                side = LEFT;
                            } else if (idx == 1) {
                                x = chan + 1;
                                y = j - 1;
                                side = TOP;
                            } else if (idx == 2) {
                                x = chan;
                                y = j;
                                side = BOTTOM;
                            } else if (idx == 3) {
                                x = chan;
                                y = j - 1;
                                side = RIGHT;
                            }
                        }
                        if(y < 0) continue;
                    }
                    /* PAJ - if the pointed to is an EMPTY then shouldn't look for ipins */
                    if (device_ctx.grid[x][y].type == device_ctx.EMPTY_TYPE)
                        continue;

                    /* Move from logical (straight) to physical (twisted) track index
                 	* - algorithm assigns ipin connections to same physical track index
                 	* so that the logical track gets distributed uniformly */

                    phy_track = vpr_to_phy_track(track, chan, j, seg_details, directionality);
                    phy_track %= tracks_per_chan;

                    /* We need the type to find the ipin map for this type */
                    type = device_ctx.grid[x][y].type;
                    int width_offset = device_ctx.grid[x][y].width_offset;
                    int height_offset = device_ctx.grid[x][y].height_offset;

                    max_conn = track_to_pin_lookup[type->index][phy_track][width_offset][height_offset][side][idx].size();
                    for (iconn = 0; iconn < max_conn; iconn++) {
                        ipin = track_to_pin_lookup[type->index][phy_track][width_offset][height_offset][side][idx][iconn];

                        /* Check there is a connection and Fc map isn't wrong */
                        /*int to_node = get_rr_node_index(L_rr_node_indices, x + width_offset, y + height_offset, IPIN, ipin, side);*/
                        int to_node = get_rr_node_index(L_rr_node_indices, x, y, IPIN, ipin, side);
                        if (to_node >= 0) {
                            rr_edges_to_create.emplace_back(from_rr_node, to_node, wire_to_ipin_switch);
                            ++num_conn;
                        }
                    }
                }
            }
        }
    }
    return (num_conn);
}

/*
 * Collects the edges fanning-out of the 'from' track which connect to the 'to'
 * tracks, according to the switch block pattern.
 *
 * It returns the number of connections added, and updates edge_list_ptr to
 * point at the head of the (extended) linked list giving the nodes to which
 * this segment connects and the switch type used to connect to each.
 *
 * An edge is added from this segment to a y-segment if:
 * (1) this segment should have a switch box at that location, or
 * (2) the y-segment to which it would connect has a switch box, and the switch
 *     type of that y-segment is unbuffered (bidirectional pass transistor).
 *
 * For bidirectional:
 * If the switch in each direction is a pass transistor (unbuffered), both
 * switches are marked as being of the types of the larger (lower R) pass
 * transistor.
 */
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
                        std::map<int, int>& bend_segment_map) {
    int to_chan, to_sb;
    std::vector<int> conn_tracks;
    bool from_is_sblock, is_behind, Fs_clipped;
    enum e_side from_side_a, from_side_b, to_side;
    bool custom_switch_block;

    /* check whether a custom switch block will be used */
    custom_switch_block = false;
    if (sb_conn_map) {
        custom_switch_block = true;
        VTR_ASSERT(switch_block_conn.empty());
    }

    VTR_ASSERT_MSG(from_seg == get_seg_start(from_seg_details, from_track, from_chan, from_seg), "From segment location must be a the wire start point");

    int from_switch = from_seg_details[from_track].arch_wire_switch();

    //The absolute coordinate along the channel where the switch block at the
    //beginning of the current wire segment is located
    int start_sb_seg = from_seg - 1;

    //The absolute coordinate along the channel where the switch block at the
    //end of the current wire segment is lcoated
    int end_sb_seg = get_seg_end(from_seg_details, from_track, from_seg, from_chan, chan_len);

    /* Figure out the sides of SB the from_wire will use */
    if (CHANX == from_type) {
        from_side_a = RIGHT;
        from_side_b = LEFT;
    } else {
        VTR_ASSERT(CHANY == from_type);
        from_side_a = TOP;
        from_side_b = BOTTOM;
    }

    //Set the loop bounds so we iterate over the whole wire segment
    int start = start_sb_seg;
    int end = end_sb_seg;

    //If source and destination segments both lie along the same channel
    //we clip the loop bounds to the switch blocks of interest and proceed
    //normally
    if (to_type == from_type) {
        start = to_seg - 1;
        end = to_seg;
    }

    //Walk along the 'from' wire segment identifying if a switchblock is located
    //at each coordinate and add any related fan-out connections to the 'from' wire segment
    int num_conn = 0;
    for (int sb_seg = start; sb_seg <= end; ++sb_seg) {
        if (sb_seg < start_sb_seg || sb_seg > end_sb_seg) {
            continue;
        }

        /* Figure out if we are at a sblock */
        from_is_sblock = is_sblock(from_chan, from_seg, sb_seg, from_track,
                                   from_seg_details, directionality);
        if (sb_seg == end_sb_seg || sb_seg == start_sb_seg) {
            /* end of wire must be an sblock */
            from_is_sblock = true;
        }

        auto switch_override = should_create_switchblock(grid, from_chan, sb_seg, from_type, to_type);
        if (switch_override == NO_SWITCH) {
            continue; //Do not create an SB here
        }

        /* Get the coordinates of the current SB from the perspective of the destination channel.
         * i.e. for segments laid in the x-direction, sb_seg corresponds to the x coordinate and from_chan to the y,
         * but for segments in the y-direction, from_chan is the x coordinate and sb_seg is the y. So here we reverse
         * the coordinates if necessary */
        if (from_type == to_type) {
            //Same channel
            to_chan = from_chan;
            to_sb = sb_seg;
        } else {
            VTR_ASSERT(from_type != to_type);
            //Different channels
            to_chan = sb_seg;
            to_sb = from_chan;
        }

        /* to_chan_details may correspond to an x-directed or y-directed channel, depending for which
         * channel type this function is used; so coordinates are reversed as necessary */
        if (to_type == CHANX) {
            to_seg_details = to_chan_details[to_seg][to_chan].data();
        } else {
            to_seg_details = to_chan_details[to_chan][to_seg].data();
        }

        if (to_seg_details[0].length() == 0)
            continue;

        /* Figure out whether the switch block at the current sb_seg coordinate is *behind*
         * the target channel segment (with respect to VPR coordinate system) */
        is_behind = false;
        if (to_type == from_type) {
            if (sb_seg == start) {
                is_behind = true;
            } else {
                is_behind = false;
            }
        } else {
            VTR_ASSERT((to_seg == from_chan) || (to_seg == (from_chan + 1)));
            if (to_seg > from_chan) {
                is_behind = true;
            }
        }

        /* Figure out which side of the SB the destination segment lies on */
        if (CHANX == to_type) {
            to_side = (is_behind ? RIGHT : LEFT);
        } else {
            VTR_ASSERT(CHANY == to_type);
            to_side = (is_behind ? TOP : BOTTOM);
        }

        /* To get to the destination seg/chan, the source track can connect to the SB from
         * one of two directions. If we're in CHANX, we can connect to it from the left or
         * right, provided we're not at a track endpoint. And similarly for a source track
         * in CHANY. */
        /* Do edges going from the right SB side (if we're in CHANX) or top (if we're in CHANY).
         * However, can't connect to right (top) if already at rightmost (topmost) track end */
        if (sb_seg < end_sb_seg) {
            if (custom_switch_block) {
                if (DEC_DIRECTION == from_seg_details[from_track].direction() || BI_DIRECTIONAL == directionality) {
                    num_conn += get_track_to_chan_seg(from_track, to_chan, to_seg,
                                                      to_type, from_side_a, to_side,
                                                      switch_override,
                                                      L_rr_node_indices,
                                                      sb_conn_map, from_rr_node, rr_edges_to_create);
                }
            } else {
                if (BI_DIRECTIONAL == directionality) {
                    /* For bidir, the target segment might have an unbuffered (bidir pass transistor)
                     * switchbox, so we follow through regardless of whether the current segment has an SB */
                    conn_tracks = switch_block_conn[from_side_a][to_side][from_track];
                    num_conn += get_bidir_track_to_chan_seg(conn_tracks,
                                                            L_rr_node_indices, to_chan, to_seg, to_sb, to_type,
                                                            to_seg_details, from_is_sblock, from_switch,
                                                            switch_override,
                                                            directionality, from_rr_node, rr_edges_to_create);
                }
                if (UNI_DIRECTIONAL == directionality) {
                    /* No fanout if no SB. */
                    /* Also, we are connecting from the top or right of SB so it
                     * makes the most sense to only get there from DEC_DIRECTION wires. */
                    if ((from_is_sblock) && (DEC_DIRECTION == from_seg_details[from_track].direction())) {
                        num_conn += get_unidir_track_to_chan_seg(from_track, to_chan,
                                                                 to_seg, to_sb, to_type, max_chan_width, grid,
                                                                 from_side_a, to_side, Fs_per_side,
                                                                 sblock_pattern,
                                                                 switch_override,
                                                                 L_rr_node_indices, to_seg_details,
                                                                 &Fs_clipped, from_rr_node, rr_edges_to_create,
                                                                 bend_delayless_switch,
                                                                 bend_segment_map);
                    }
                }
            }
        }

        /* Do the edges going from the left SB side (if we're in CHANX) or bottom (if we're in CHANY)
         * However, can't connect to left (bottom) if already at leftmost (bottommost) track end */
        if (sb_seg > start_sb_seg) {
            if (custom_switch_block) {
                if (INC_DIRECTION == from_seg_details[from_track].direction() || BI_DIRECTIONAL == directionality) {
                    num_conn += get_track_to_chan_seg(from_track, to_chan, to_seg,
                                                      to_type, from_side_b, to_side,
                                                      switch_override,
                                                      L_rr_node_indices,
                                                      sb_conn_map, from_rr_node, rr_edges_to_create);
                }
            } else {
                if (BI_DIRECTIONAL == directionality) {
                    /* For bidir, the target segment might have an unbuffered (bidir pass transistor)
                     * switchbox, so we follow through regardless of whether the current segment has an SB */
                    conn_tracks = switch_block_conn[from_side_b][to_side][from_track];
                    num_conn += get_bidir_track_to_chan_seg(conn_tracks,
                                                            L_rr_node_indices, to_chan, to_seg, to_sb, to_type,
                                                            to_seg_details, from_is_sblock, from_switch,
                                                            switch_override,
                                                            directionality, from_rr_node, rr_edges_to_create);
                }
                if (UNI_DIRECTIONAL == directionality) {
                    /* No fanout if no SB. */
                    /* Also, we are connecting from the bottom or left of SB so it
                     * makes the most sense to only get there from INC_DIRECTION wires. */
                    if ((from_is_sblock)
                        && (INC_DIRECTION == from_seg_details[from_track].direction())) {
                        num_conn += get_unidir_track_to_chan_seg(from_track, to_chan,
                                                                 to_seg, to_sb, to_type, max_chan_width, grid,
                                                                 from_side_b, to_side, Fs_per_side,
                                                                 sblock_pattern,
                                                                 switch_override,
                                                                 L_rr_node_indices, to_seg_details,
                                                                 &Fs_clipped, from_rr_node, rr_edges_to_create,
                                                                 bend_delayless_switch,
                                                                 bend_segment_map);
                    }
                }
            }
        }
    }

    return num_conn;
}

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
                                        std::map<int, int>& bend_segment_map) {
    int to_chan, to_sb;
    std::vector<int> conn_tracks;
    bool is_behind, Fs_clipped;
    enum e_side from_side_a, from_side_b, to_side;

    //The absolute coordinate along the channel where the switch block at the
    //beginning of the current wire segment is located
    int start_sb_seg = from_seg - 1;

    //The absolute coordinate along the channel where the switch block at the
    //end of the current wire segment is lcoated
    int end_sb_seg = get_seg_end(from_seg_details, from_track, from_seg, from_chan, chan_len);

    /* Figure out the sides of SB the from_wire will use */
    if (CHANX == from_type) {
        from_side_a = RIGHT;
        from_side_b = LEFT;
    } else {
        VTR_ASSERT(CHANY == from_type);
        from_side_a = TOP;
        from_side_b = BOTTOM;
    }

    //Walk along the 'from' wire segment identifying if a switchblock is located
    //at each coordinate and add any related fan-out connections to the 'from' wire segment
    int num_conn = 0;
    
    if (sb_seg < start_sb_seg || sb_seg > end_sb_seg) return num_conn;
    
    auto switch_override = should_create_switchblock(grid, from_chan, sb_seg, from_type, to_type);
    if (switch_override == NO_SWITCH) {
        return num_conn; //Do not create an SB here
    }

    /* Get the coordinates of the current SB from the perspective of the destination channel.
     * i.e. for segments laid in the x-direction, sb_seg corresponds to the x coordinate and from_chan to the y,
     * but for segments in the y-direction, from_chan is the x coordinate and sb_seg is the y. So here we reverse
     * the coordinates if necessary */
    if (from_type == to_type) {
        //Same channel
        to_chan = from_chan;
        to_sb = sb_seg;
    } else {
        VTR_ASSERT(from_type != to_type);
        //Different channels
        to_chan = sb_seg;
        to_sb = from_chan;
    }

    /* to_chan_details may correspond to an x-directed or y-directed channel, depending for which
     * channel type this function is used; so coordinates are reversed as necessary */
    if (to_type == CHANX) {
        to_seg_details = to_chan_details[to_seg][to_chan].data();
    } else {
        to_seg_details = to_chan_details[to_chan][to_seg].data();
    }

    if (to_seg_details[0].length() == 0)
        return num_conn;

    /* Figure out whether the switch block at the current sb_seg coordinate is *behind*
     * the target channel segment (with respect to VPR coordinate system) */
    is_behind = false;
    if (to_type == from_type) {
        if (sb_seg == to_seg - 1) {
            is_behind = true;
        } else {
            is_behind = false;
        }
    } else {
        VTR_ASSERT((to_seg == from_chan) || (to_seg == (from_chan + 1)));
        if (to_seg > from_chan) {
            is_behind = true;
        }
    }

    /* Figure out which side of the SB the destination segment lies on */
    if (CHANX == to_type) {
        to_side = (is_behind ? RIGHT : LEFT);
    } else {
        VTR_ASSERT(CHANY == to_type);
        to_side = (is_behind ? TOP : BOTTOM);
    }

    /* To get to the destination seg/chan, the source track can connect to the SB from
     * one of two directions. If we're in CHANX, we can connect to it from the left or
     * right, provided we're not at a track endpoint. And similarly for a source track
     * in CHANY. */
    /* Do edges going from the right SB side (if we're in CHANX) or top (if we're in CHANY).
     * However, can't connect to right (top) if already at rightmost (topmost) track end */
    if (sb_seg < end_sb_seg) {
        if (UNI_DIRECTIONAL == directionality) {
            /* No fanout if no SB. */
            /* Also, we are connecting from the top or right of SB so it
              * makes the most sense to only get there from DEC_DIRECTION wires. */
            if (DEC_DIRECTION == from_seg_details[from_track].direction()) {
                num_conn += get_unidir_track_to_chan_seg(from_track, to_chan,
                                                             to_seg, to_sb, to_type, max_chan_width, grid,
                                                             from_side_a, to_side, Fs_per_side,
                                                             sblock_pattern,
                                                             switch_override,
                                                             L_rr_node_indices, to_seg_details,
                                                             &Fs_clipped, from_rr_node, rr_edges_to_create,
                                                             bend_delayless_switch, bend_segment_map);
            }
        }
        else{
            VPR_THROW(VPR_ERROR_ROUTE,
                      "in get_track_to_tracks_for_gsb: only support UNI_DIRECTIONAL wire!");
        }
    }

    /* Do the edges going from the left SB side (if we're in CHANX) or bottom (if we're in CHANY)
     * However, can't connect to left (bottom) if already at leftmost (bottommost) track end */
    if (sb_seg > start_sb_seg) {
        if (UNI_DIRECTIONAL == directionality) {
            /* No fanout if no SB. */
            /* Also, we are connecting from the bottom or left of SB so it
             * makes the most sense to only get there from INC_DIRECTION wires. */
            if (INC_DIRECTION == from_seg_details[from_track].direction()) {
                num_conn += get_unidir_track_to_chan_seg(from_track, to_chan,
                                                         to_seg, to_sb, to_type, max_chan_width, grid,
                                                         from_side_b, to_side, Fs_per_side,
                                                         sblock_pattern,
                                                         switch_override,
                                                         L_rr_node_indices, to_seg_details,
                                                         &Fs_clipped, from_rr_node, rr_edges_to_create,
                                                         bend_delayless_switch, bend_segment_map);
            }
        }
        else{
            VPR_THROW(VPR_ERROR_ROUTE,
                      "in get_track_to_tracks_for_gsb: only support UNI_DIRECTIONAL wire!");
        }
    }
    return num_conn;
}

static int get_bidir_track_to_chan_seg(const std::vector<int> conn_tracks,
                                       const t_rr_node_indices& L_rr_node_indices,
                                       const int to_chan,
                                       const int to_seg,
                                       const int to_sb,
                                       const t_rr_type to_type,
                                       const t_chan_seg_details* seg_details,
                                       const bool from_is_sblock,
                                       const int from_switch,
                                       const int switch_override,
                                       const enum e_directionality directionality,
                                       const int from_rr_node,
                                       t_rr_edge_info_set& rr_edges_to_create) {
    unsigned iconn;
    int to_track, to_node, to_switch, num_conn, to_x, to_y, i;
    bool to_is_sblock;
    short switch_types[2];

    /* x, y coords for get_rr_node lookups */
    if (CHANX == to_type) {
        to_x = to_seg;
        to_y = to_chan;
    } else {
        VTR_ASSERT(CHANY == to_type);
        to_x = to_chan;
        to_y = to_seg;
    }

    /* Go through the list of tracks we can connect to */
    num_conn = 0;
    for (iconn = 0; iconn < conn_tracks.size(); ++iconn) {
        to_track = conn_tracks[iconn];
        to_node = get_rr_node_index(L_rr_node_indices, to_x, to_y, to_type, to_track);

        if (to_node == OPEN) {
            continue;
        }

        /* Get the switches for any edges between the two tracks */
        to_switch = seg_details[to_track].arch_wire_switch();

        to_is_sblock = is_sblock(to_chan, to_seg, to_sb, to_track, seg_details,
                                 directionality);
        get_switch_type(from_is_sblock, to_is_sblock, from_switch, to_switch,
                        switch_override,
                        switch_types);

        /* There are up to two switch edges allowed from track to track */
        for (i = 0; i < 2; ++i) {
            /* If the switch_type entry is empty, skip it */
            if (OPEN == switch_types[i]) {
                continue;
            }

            /* Add the edge to the list */
            rr_edges_to_create.emplace_back(from_rr_node, to_node, switch_types[i]);
            ++num_conn;
        }
    }

    return num_conn;
}

/* Figures out the edges that should connect the given wire segment to the given
 * channel segment, adds these edges to 'edge_list' and returns the number of
 * edges added .
 * See route/build_switchblocks.c for a detailed description of how the switch block
 * connection map sb_conn_map is generated. */
static int get_track_to_chan_seg(const int from_wire,
                                 const int to_chan,
                                 const int to_seg,
                                 const t_rr_type to_chan_type,
                                 const e_side from_side,
                                 const e_side to_side,
                                 const int switch_override,
                                 const t_rr_node_indices& L_rr_node_indices,
                                 t_sb_connection_map* sb_conn_map,
                                 const int from_rr_node,
                                 t_rr_edge_info_set& rr_edges_to_create) {
    int edge_count = 0;
    int to_x, to_y;
    int tile_x, tile_y;

    /* get x/y coordinates from seg/chan coordinates */
    if (CHANX == to_chan_type) {
        to_x = tile_x = to_seg;
        to_y = tile_y = to_chan;
        if (RIGHT == to_side) {
            tile_x--;
        }
    } else {
        VTR_ASSERT(CHANY == to_chan_type);
        to_x = tile_x = to_chan;
        to_y = tile_y = to_seg;
        if (TOP == to_side) {
            tile_y--;
        }
    }

    /* get coordinate to index into the SB map */
    Switchblock_Lookup sb_coord(tile_x, tile_y, from_side, to_side);
    if (sb_conn_map->count(sb_coord) > 0) {
        /* get reference to the connections vector which lists all destination wires for a given source wire
         * at a specific coordinate sb_coord */
        vector<t_switchblock_edge>& conn_vector = (*sb_conn_map)[sb_coord];

        /* go through the connections... */
        for (int iconn = 0; iconn < (int)conn_vector.size(); ++iconn) {
            if (conn_vector.at(iconn).from_wire != from_wire) continue;

            int to_wire = conn_vector.at(iconn).to_wire;
            int to_node = get_rr_node_index(L_rr_node_indices, to_x, to_y, to_chan_type, to_wire);

            if (to_node == OPEN) {
                continue;
            }

            /* Get the index of the switch connecting the two wires */
            int src_switch = conn_vector[iconn].switch_ind;

            //Apply any switch overrides
            if (should_apply_switch_override(switch_override)) {
                src_switch = switch_override;
            }

            rr_edges_to_create.emplace_back(from_rr_node, to_node, src_switch);
            ++edge_count;

            auto& device_ctx = g_vpr_ctx.device();

            if (device_ctx.arch_switch_inf[src_switch].directionality() == BI_DIRECTIONAL) {
                //Add reverse edge since bi-directional
                rr_edges_to_create.emplace_back(to_node, from_rr_node, src_switch);
                ++edge_count;
            }
        }
    } else {
        /* specified sb_conn_map entry does not exist -- do nothing */
    }
    return edge_count;
}

static int get_unidir_track_to_chan_seg(const int from_track,
                                        const int to_chan,
                                        const int to_seg,
                                        const int to_sb,
                                        const t_rr_type to_type,
                                        const int max_chan_width,
                                        const DeviceGrid& grid,
                                        const enum e_side from_side,
                                        const enum e_side to_side,
                                        const int Fs_per_side,
                                        t_sblock_pattern& sblock_pattern,
                                        const int switch_override,
                                        const t_rr_node_indices& L_rr_node_indices,
                                        const t_chan_seg_details* seg_details,
                                        bool* Fs_clipped,
                                        const int from_rr_node,
                                        t_rr_edge_info_set& rr_edges_to_create,
                                        const int bend_delayless_switch,
                                        map<int, int>& bend_segment_map) {
    int num_labels = 0;
    int* mux_labels = nullptr;

    /* x, y coords for get_rr_node lookups */
    int to_x = (CHANX == to_type ? to_seg : to_chan);
    int to_y = (CHANX == to_type ? to_chan : to_seg);
    int sb_x = (CHANX == to_type ? to_sb : to_chan);
    int sb_y = (CHANX == to_type ? to_chan : to_sb);
    int max_len = (CHANX == to_type ? grid.width() : grid.height()) - 2; //-2 for no perimeter channels

    enum e_direction to_dir = DEC_DIRECTION;
    if (to_sb < to_seg) {
        to_dir = INC_DIRECTION;
    }

    *Fs_clipped = false;

    /* get list of muxes to which we can connect */
    int dummy;
    bool is_dangling_side = label_dangling_side(to_side, sb_x, sb_y, grid);
    mux_labels = label_wire_muxes(to_chan, to_seg, seg_details, UNDEFINED, max_len,
                                  to_dir, max_chan_width, false, &num_labels, &dummy, is_dangling_side);

    /* Can't connect if no muxes. */
    if (num_labels < 1) {
        if (mux_labels) {
            vtr::free(mux_labels);
            mux_labels = nullptr;
        }
        return 0;
    }

    /* Check if Fs demand was too high. */
    if (Fs_per_side > num_labels) {
        *Fs_clipped = true;
    }

    /* Handle Fs > 3 by assigning consecutive muxes. */
    int count = 0;
    for (int i = 0; i < Fs_per_side; ++i) {
        /* Get the target label */
        for (int j = 0; j < 4; j = j + 2) {
            /* Use the balanced labeling for passing and fringe wires */
            int to_mux = sblock_pattern[sb_x][sb_y][from_side][to_side][from_track][j];
            int to_track = sblock_pattern[sb_x][sb_y][from_side][to_side][from_track][j + 1];
            if (to_mux == UN_SET && to_track == UN_SET)
                continue;

            if (to_mux == UN_SET && to_track != UN_SET) {
                /* Add edge to list. */

                int to_node = get_rr_node_index(L_rr_node_indices, to_x, to_y, to_type, to_track);

                if (to_node == OPEN) {
                    continue;
                }

                /*if(bend_segment_map.count(to_node) != 0){
                    int a=0;
                    a++;
                }*/
                VTR_ASSERT(bend_segment_map.count(to_node) == 0);
                bend_segment_map[to_node] = from_rr_node;
                rr_edges_to_create.emplace_back(from_rr_node, to_node, bend_delayless_switch);
                ++count;
                continue;
            }

            //int to_track = sblock_pattern[sb_x][sb_y][from_side][to_side][from_track][j + 1];
            if (to_track == UN_SET) {
                to_track = mux_labels[(to_mux + i) % num_labels];
                sblock_pattern[sb_x][sb_y][from_side][to_side][from_track][j + 1] = to_track;
            }

            int to_node = get_rr_node_index(L_rr_node_indices, to_x, to_y, to_type, to_track);

            if (to_node == OPEN) {
                continue;
            }

            // check if the middle bend segment should have a mux
            if (seg_details[to_track].part_idx() > 0)
                VTR_ASSERT(is_dangling_side);

            //Determine which switch to use
            int iswitch = seg_details[to_track].arch_wire_switch();
            if (iswitch == OPEN && seg_details[to_track].part_idx() > 0)
                iswitch = seg_details[to_track].first_wire_switch();

            //Apply any switch overrides
            if (should_apply_switch_override(switch_override)) {
                iswitch = switch_override;
            }
            VTR_ASSERT(iswitch != OPEN);

            /* Add edge to list. */
            rr_edges_to_create.emplace_back(from_rr_node, to_node, iswitch);
            ++count;

            auto& device_ctx = g_vpr_ctx.device();
            if (device_ctx.arch_switch_inf[iswitch].directionality() == BI_DIRECTIONAL) {
                //Add reverse edge since bi-directional
                rr_edges_to_create.emplace_back(to_node, from_rr_node, iswitch);
                ++count;
            }
        }
    }

    if (mux_labels) {
        vtr::free(mux_labels);
        mux_labels = nullptr;
    }
    return count;
}

bool is_sblock(const int chan, int wire_seg, const int sb_seg, const int track, const t_chan_seg_details* seg_details, const enum e_directionality directionality) {
    int length, ofs, fac;

    fac = 1;
    if (UNI_DIRECTIONAL == directionality) {
        fac = 2;
    }

    length = seg_details[track].length();

    /* Make sure they gave us correct start */
    wire_seg = get_seg_start(seg_details, track, chan, wire_seg);

    ofs = sb_seg - wire_seg + 1; /* Ofset 0 is behind us, so add 1 */

    VTR_ASSERT(ofs >= 0);
    VTR_ASSERT(ofs < (length + 1));

    /* If unidir segment that is going backwards, we need to flip the ofs */
    if ((track % fac) > 0) {
        ofs = length - ofs;
    }

    return seg_details[track].sb(ofs);
}

static void get_switch_type(bool is_from_sblock,
                            bool is_to_sblock,
                            short from_node_switch,
                            short to_node_switch,
                            const int switch_override,
                            short switch_types[2]) {
    /* This routine looks at whether the from_node and to_node want a switch,  *
     * and what type of switch is used to connect *to* each type of node       *
     * (from_node_switch and to_node_switch).  It decides what type of switch, *
     * if any, should be used to go from from_node to to_node.  If no switch   *
     * should be inserted (i.e. no connection), it returns OPEN.  Its returned *
     * values are in the switch_types array.  It needs to return an array      *
     * because some topologies (e.g. bi-dir pass gates) result in two switches.*/

    auto& device_ctx = g_vpr_ctx.device();

    switch_types[0] = NO_SWITCH;
    switch_types[1] = NO_SWITCH;

    if (switch_override == NO_SWITCH) {
        return; //No switches
    }

    if (should_apply_switch_override(switch_override)) {
        //Use the override switches instead
        from_node_switch = switch_override;
        to_node_switch = switch_override;
    }

    int used = 0;
    bool forward_switch = false;
    bool backward_switch = false;

    /* Connect forward if we are a sblock */
    if (is_from_sblock) {
        switch_types[used] = to_node_switch;
        ++used;

        forward_switch = true;
    }

    /* Check for reverse switch */
    if (is_to_sblock) {
        if (device_ctx.arch_switch_inf[from_node_switch].directionality() == e_directionality::BI_DIRECTIONAL) {
            switch_types[used] = from_node_switch;
            ++used;

            backward_switch = true;
        }
    }

    /* Take the larger switch if there are two of the same type */
    if (forward_switch
        && backward_switch
        && (device_ctx.arch_switch_inf[from_node_switch].type() == device_ctx.arch_switch_inf[to_node_switch].type())) {
        //Sanity checks
        VTR_ASSERT_SAFE_MSG(device_ctx.arch_switch_inf[from_node_switch].type() == device_ctx.arch_switch_inf[to_node_switch].type(), "Same switch type");
        VTR_ASSERT_MSG(device_ctx.arch_switch_inf[to_node_switch].directionality() == e_directionality::BI_DIRECTIONAL, "Bi-dir to switch");
        VTR_ASSERT_MSG(device_ctx.arch_switch_inf[from_node_switch].directionality() == e_directionality::BI_DIRECTIONAL, "Bi-dir from switch");

        /* Take the smaller index unless the other
         * switch is bigger (smaller R). */

        int first_switch = min(to_node_switch, from_node_switch);
        int second_switch = max(to_node_switch, from_node_switch);

        if (used < 2) {
            VPR_THROW(VPR_ERROR_ROUTE,
                      "Expected 2 switches (forward and back) between RR nodes (found %d switches, min switch index: %d max switch index: %d)",
                      used, first_switch, second_switch);
        }

        int switch_to_use = first_switch;
        if (device_ctx.arch_switch_inf[second_switch].R < device_ctx.arch_switch_inf[first_switch].R) {
            switch_to_use = second_switch;
        }

        for (int i = 0; i < used; ++i) {
            switch_types[i] = switch_to_use;
        }
    }
}

static int vpr_to_phy_track(const int itrack,
                            const int chan_num,
                            const int seg_num,
                            const t_chan_seg_details* seg_details,
                            const enum e_directionality directionality) {
    int group_start, group_size;
    int vpr_offset_for_first_phy_track;
    int vpr_offset, phy_offset;
    int phy_track;
    int fac;

    /* Assign in pairs if unidir. */
    fac = 1;
    if (UNI_DIRECTIONAL == directionality) {
        fac = 2;
    }

    group_start = seg_details[itrack].group_start();
    group_size = seg_details[itrack].group_size();

    vpr_offset_for_first_phy_track = (chan_num + seg_num - 1)
                                     % (group_size / fac);
    vpr_offset = (itrack - group_start) / fac;
    phy_offset = (vpr_offset_for_first_phy_track + vpr_offset)
                 % (group_size / fac);
    phy_track = group_start + (fac * phy_offset) + (itrack - group_start) % fac;

    return phy_track;
}

t_sblock_pattern alloc_sblock_pattern_lookup(const DeviceGrid& grid,
                                             const int max_chan_width) {
    /* loading up the sblock connection pattern matrix. It's a huge matrix because
     * for nonquantized W, it's impossible to make simple permutations to figure out
     * where muxes are and how to connect to them such that their sizes are balanced */

    /* Do chunked allocations to make freeing easier, speed up malloc and free, and
     * reduce some of the memory overhead. Could use fewer malloc's but this way
     * avoids all considerations of pointer sizes and allignment. */

    /* Alloc each list of pointers in one go. items is a running product that increases
     * with each new dimension of the matrix. */

    VTR_ASSERT(grid.width() > 0);
    VTR_ASSERT(grid.height() > 0);
    VTR_ASSERT(max_chan_width >= 0);

    t_sblock_pattern sblock_pattern({{
                                        grid.width() - 1,
                                        grid.height() - 1,
                                        4, //From side
                                        4, //To side
                                        size_t(max_chan_width),
                                        4 //to_mux, to_trac, alt_mux, alt_track
                                    }},
                                    UN_SET);

    /* This is the outer pointer to the full matrix */
    return sblock_pattern;
}

static short find_to_side_bend_point(int itrack, const t_chan_seg_details* from_seg_details, const t_chan_seg_details* to_seg_details, const int* incoming_end_label, const int num_incoming_bend_end, const int* wire_bend_label, const int num_bend_wires) {
    std::vector<int> from_points;
    std::vector<int> to_points;
    int i;
    short to_track;

    VTR_ASSERT(num_bend_wires == num_incoming_bend_end);
    for (i = 0; i < num_bend_wires; i++) {
        int wire_idx = wire_bend_label[i];
        if (to_seg_details[wire_idx].type_name() == from_seg_details[itrack].type_name()) {
            if (to_seg_details[wire_idx].part_idx() == from_seg_details[itrack].part_idx() + 1)
                to_points.push_back(wire_idx);
        }
    }

    for (i = 0; i < num_incoming_bend_end; i++) {
        int wire_idx = incoming_end_label[i];
        if (from_seg_details[wire_idx].type_name() == from_seg_details[itrack].type_name()) {
            if (from_seg_details[wire_idx].part_idx() == from_seg_details[itrack].part_idx())
                from_points.push_back(wire_idx);
        }
    }

    VTR_ASSERT(!from_points.empty());
    VTR_ASSERT(!to_points.empty());
    VTR_ASSERT(from_points.size() == to_points.size());

    auto itr = find(from_points.begin(), from_points.end(), itrack);
    /*if(itr == from_points.end()){
        int a = 0;
        a++;
    }*/
    VTR_ASSERT(itr != from_points.end());
    auto idx_from = std::distance(std::begin(from_points), itr);
    to_track = (short)to_points[idx_from];

    return to_track;
}

void load_sblock_bend_pattern_lookup_for_gsb(const int i,
                                             const int j,
                                             const DeviceGrid& grid,
                                             const t_chan_width* nodes_per_chan,
                                             const t_chan_details& chan_details_x,
                                             const t_chan_details& chan_details_y,
                                             const int /*Fs*/,
                                             const enum e_switch_block_type switch_block_type,
                                             t_sblock_pattern& sblock_pattern) {
    /* This routine loads a lookup table for sblock topology. The lookup table is huge
     * because the sblock varies from location to location. The i, j means the owning
     * location of the sblock under investigation. */

    /* SB's have coords from (0, 0) to (grid.width()-2, grid.height()-2) */
    VTR_ASSERT(i >= 0);
    VTR_ASSERT(i <= int(grid.width()) - 2);
    VTR_ASSERT(j >= 0);
    VTR_ASSERT(j <= int(grid.height()) - 2);

    /* May 12 - 15, 2007
     *
     * I identify three types of sblocks in the chip: 1) The core sblock, whose special
     * property is that the number of muxes (and ending wires) on each side is the same (very useful
     * property, since it leads to a N-to-N assignment problem with ending wires). 2) The corner sblock
     * which is same as a L=1 core sblock with 2 sides only (again N-to-N assignment problem). 3) The
     * fringe / chip edge sblock which is most troublesome, as balance in each side of muxes is
     * attainable but balance in the entire sblock is not. The following code first identifies the
     * incoming wires, which can be classified into incoming passing wires with sblock and incoming
     * ending wires (the word "incoming" is sometimes dropped for ease of discussion). It appropriately
     * labels all the wires on each side by the following order: By the call to label_incoming_wires,
     * which labels for one side, the order is such that the incoming ending wires (always with sblock)
     * are labelled first 0,1,2,... p-1, then the incoming passing wires with sblock are labelled
     * p,p+1,p+2,... k-1 (for total of k). By this convention, one can easily distinguish the ending
     * wires from the passing wires by checking a label against num_ending_wires variable.
     *
     * After labelling all the incoming wires, this routine labels the muxes on the side we're currently
     * connecting to (iterated for four sides of the sblock), called the to_side. The label scheme is
     * the natural order of the muxes by their track #. Also we find the number of muxes.
     *
     * For each to_side, the total incoming wires that connect to the muxes on to_side
     * come from three sides: side_1 (to_side's right), side_2 (to_side's left) and opp_side.
     * The problem of balancing mux size is then: considering all incoming passing wires
     * with sblock on side_1, side_2 and opp_side, how to assign them to the muxes on to_side
     * (with specified Fs) in a way that mux size is imbalanced by at most 1. I solve this
     * problem by this approach: the first incoming passing wire will connect to 0, 1, 2,
     * ..., Fs_per_side - 1, then the next incoming passing wire will connect to
     * Fs_per_side, Fs_per_side+1, ..., Fs_per_side*2-1, and so on. This consistent STAGGERING
     * ensures N-to-N assignment is perfectly balanced and M-to-N assignment is imbalanced by no
     * more than 1.
     */

    /* SB's range from (0, 0) to (grid.width() - 2, grid.height() - 2) */
    /* First find all four sides' incoming wires */
    int* wire_mux_on_track[4];
    int* incoming_wire_label[4];
    int num_incoming_wires[4];
    int num_ending_wires[4];
    int num_wire_muxes[4];

    int* wire_bend_label[4];
    int num_bend_wires[4];
    int* incoming_bend_end_label[4];
    int num_incoming_bend[4];

    bool is_dangling_side[4];

    /* "Label" if the side of switch block is the dangling side according to the location of switch block */
    for (e_side side : {TOP, RIGHT, BOTTOM, LEFT}) {
        is_dangling_side[side] = label_dangling_side(side, i, j, grid);
    }

    /* "Label" the wires around the switch block by connectivity. */
    for (e_side side : {TOP, RIGHT, BOTTOM, LEFT}) {
        /* Assume the channel segment doesn't exist. */
        wire_mux_on_track[side] = nullptr;
        incoming_wire_label[side] = nullptr;
        num_incoming_wires[side] = 0;
        num_ending_wires[side] = 0;
        num_wire_muxes[side] = 0;

        wire_bend_label[side] = nullptr;
        num_bend_wires[side] = 0;
        incoming_bend_end_label[side] = nullptr;
        num_incoming_bend[side] = 0;

        /* Skip the side and leave the zero'd value if the
         * channel segment doesn't exist. */
        bool skip = true;
        switch (side) {
            case TOP:
                if (j < int(grid.height()) - 2) {
                    skip = false;
                }
                break;
            case RIGHT:
                if (i < int(grid.width()) - 2) {
                    skip = false;
                }
                break;
            case BOTTOM:
                if (j > 0) {
                    skip = false;
                }
                break;
            case LEFT:
                if (i > 0) {
                    skip = false;
                }
                break;
            default:
                VTR_ASSERT_MSG(false, "Unrecognzied side");
                break;
        }
        if (skip) {
            continue;
        }

        /* Figure out the channel and segment for a certain direction */
        bool vert = ((side == TOP) || (side == BOTTOM));
        bool pos_dir = ((side == TOP) || (side == RIGHT));
        int chan_len = (vert ? grid.height() : grid.width()) - 2; //-2 for no perim channels
        int chan = (vert ? i : j);
        int sb_seg = (vert ? j : i);
        int seg = (pos_dir ? (sb_seg + 1) : sb_seg);

        const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg] : chan_details_x[seg][chan]).data();
        if (seg_details[0].length() <= 0)
            continue;

        /* Figure out all the tracks on a side that are ending and the
         * ones that are passing through and have a SB. */
        enum e_direction end_dir = (pos_dir ? DEC_DIRECTION : INC_DIRECTION);
        
        incoming_wire_label[side] = label_incoming_wires(chan, seg, sb_seg,
                                                         seg_details, chan_len, end_dir, nodes_per_chan->max,
                                                         &num_incoming_wires[side], &num_ending_wires[side]);
        
        incoming_bend_end_label[side] = label_incoming_bend_wires_ending(chan, seg,
                                                                         seg_details, UNDEFINED, chan_len, end_dir, nodes_per_chan->max,
                                                                         &num_incoming_bend[side], is_dangling_side[side]);

        /* Figure out all the tracks on a side that are starting. */
        int dummy;
        enum e_direction start_dir = (pos_dir ? INC_DIRECTION : DEC_DIRECTION);
        wire_mux_on_track[side] = label_wire_muxes(chan, seg,
                                                   seg_details, UNDEFINED, chan_len, start_dir, nodes_per_chan->max,
                                                   false, &num_wire_muxes[side], &dummy, is_dangling_side[side]);

        wire_bend_label[side] = label_wire_bends(chan, seg,
                                                 seg_details, UNDEFINED, chan_len, start_dir, nodes_per_chan->max,
                                                 &num_bend_wires[side], is_dangling_side[side]);
    }

    for (e_side to_side : {TOP, RIGHT, BOTTOM, LEFT}) {
        /* Can't do anything if no muxes on this side. */
        if (num_wire_muxes[to_side] == 0)
            continue;

        /* Figure out side rotations */
        VTR_ASSERT((TOP == 0) && (RIGHT == 1) && (BOTTOM == 2) && (LEFT == 3));
        int side_cw = (to_side + 1) % 4;
        int side_opp = (to_side + 2) % 4;
        int side_ccw = (to_side + 3) % 4;

        /* Figure out the channel and segment for a certain direction */
        bool vert = ((to_side == TOP) || (to_side == BOTTOM));
        bool pos_dir = ((to_side == TOP) || (to_side == RIGHT));
        int chan = (vert ? i : j);
        int sb_seg = (vert ? j : i);
        int seg = (pos_dir ? (sb_seg + 1) : sb_seg);
        int chan_len = (vert ? grid.height() : grid.width()) - 2; //-2 for no perim channels

        const t_chan_seg_details* to_seg_details = (vert ? chan_details_y[chan][seg].data() : chan_details_x[seg][chan].data());
        if (to_seg_details[0].length() <= 0)
            continue;

        /* For the core sblock:
         * The new order for passing wires should appear as
         * 0,1,2..,scw-1, for passing wires with sblock on side_cw
         * scw,scw+1,...,sccw-1, for passing wires with sblock on side_ccw
         * sccw,sccw+1,... for passing wires with sblock on side_opp.
         * This way, I can keep the imbalance to at most 1.
         *
         * For the fringe sblocks, I don't distinguish between
         * passing and ending wires so the above statement still holds
         * if you replace "passing" by "incoming" */

        if (incoming_wire_label[side_cw]) {
            /* Figure out the channel and segment for a certain direction */
            vert = ((side_cw == TOP) || (side_cw == BOTTOM));
            pos_dir = ((side_cw == TOP) || (side_cw == RIGHT));
            chan_len = (vert ? grid.height() : grid.width()) - 2; //-2 for no perim channels
            chan = (vert ? i : j);
            sb_seg = (vert ? j : i);
            seg = (pos_dir ? (sb_seg + 1) : sb_seg);

            const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg].data() : chan_details_x[seg][chan].data());
            if (seg_details[0].length() <= 0)
                continue;

            for (int ichan = 0; ichan < nodes_per_chan->max; ichan++) { //nodes_per_chan->max == max_chan_width 不含中间节点
                int itrack = ichan;
                if (side_cw == TOP || side_cw == BOTTOM) {
                    itrack = ichan % nodes_per_chan->y_list[i];
                } else if (side_cw == RIGHT || side_cw == LEFT) {
                    itrack = ichan % nodes_per_chan->x_list[j];
                }

                /* DOWN stype from side_cw , if incoming wire is the ending of bend segment,
                 * and the bend style is the DOWN_STYLE
                 * then finding the track map bewtween side_cw and to_side
                 */

                // bend ending point
                e_direction check_dir = (pos_dir ? DEC_DIRECTION : INC_DIRECTION);
                if (incoming_wire_label[side_cw][itrack] != UN_SET
                    && incoming_wire_label[side_cw][itrack] < num_ending_wires[side_cw]
                    && seg_details[itrack].direction() == check_dir) {
                    // not the dangling incoming segment
                    int true_seg_start = seg - (seg + seg_details[itrack].length() + chan - seg_details[itrack].start()) % seg_details[itrack].length();
                    int true_seg_end = true_seg_start + seg_details[itrack].length() - 1;
                    int start = get_seg_start(seg_details, itrack, chan, seg);
                    int end = get_seg_end(seg_details, itrack, start, chan, chan_len);
                    bool clipped = (check_dir == DEC_DIRECTION ? true_seg_start != start : true_seg_end != end);
                    if ((!is_dangling_side[side_cw]) || !clipped)
                        // middle part of the bend segment
                        if (seg_details[itrack].part_idx() != seg_details[itrack].bend_len() - 1) {
                            // if bend style is the DOWN_TYPE
                            if (seg_details[itrack].switch_dir_type() == DOWN_TYPE) {
                                // assign to_track
                                sblock_pattern[i][j][side_cw][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                            incoming_bend_end_label[side_cw], num_incoming_bend[side_cw],
                                                                                                            wire_bend_label[to_side], num_bend_wires[to_side]);
                            }
                            else if (seg_details[itrack].switch_dir_type() == UP_TYPE && is_corner(grid, i, j)) {
                                sblock_pattern[i][j][side_cw][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                            incoming_bend_end_label[side_cw], num_incoming_bend[side_cw],
                                                                                                            wire_bend_label[to_side], num_bend_wires[to_side]);
                            }
                        }
                }
            }
        }

        if (incoming_wire_label[side_ccw]) {
            /* Figure out the channel and segment for a certain direction */
            vert = ((side_ccw == TOP) || (side_ccw == BOTTOM));
            pos_dir = ((side_ccw == TOP) || (side_ccw == RIGHT));
            chan = (vert ? i : j);
            sb_seg = (vert ? j : i);
            seg = (pos_dir ? (sb_seg + 1) : sb_seg);

            const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg].data() : chan_details_x[seg][chan].data());
            if (seg_details[0].length() <= 0)
                continue;

            for (int ichan = 0; ichan < nodes_per_chan->max; ichan++) {
                int itrack = ichan;
                if (side_ccw == TOP || side_ccw == BOTTOM) {
                    itrack = ichan % nodes_per_chan->y_list[i];
                } else if (side_ccw == RIGHT || side_ccw == LEFT) {
                    itrack = ichan % nodes_per_chan->x_list[j];
                }

                /* UP stype from side_ccw , if incoming wire is the ending of bend segment,
                 * and the bend style is the UP_STYPE
                 * then finding the track map bewtween side_ccw and to_side
                 */

                // bend ending point
                e_direction check_dir = (pos_dir ? DEC_DIRECTION : INC_DIRECTION);
                if (incoming_wire_label[side_ccw][itrack] != UN_SET
                    && incoming_wire_label[side_ccw][itrack] < num_ending_wires[side_ccw]
                    && seg_details[itrack].direction() == check_dir) {
                    // not the dangling incoming segment
                    int true_seg_start = seg - (seg + seg_details[itrack].length() + chan - seg_details[itrack].start()) % seg_details[itrack].length();
                    int true_seg_end = true_seg_start + seg_details[itrack].length() - 1;
                    int start = get_seg_start(seg_details, itrack, chan, seg);
                    int end = get_seg_end(seg_details, itrack, start, chan, chan_len);
                    bool clipped = (check_dir == DEC_DIRECTION ? true_seg_start != start : true_seg_end != end);
                    if (!clipped || (!is_dangling_side[side_ccw]))
                        // middle part of the bend segment
                        if (seg_details[itrack].part_idx() != seg_details[itrack].bend_len() - 1) {
                            // if bend style is the UP_TYPE
                            if (seg_details[itrack].switch_dir_type() == UP_TYPE) {
                                // assign to_track
                                sblock_pattern[i][j][side_ccw][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                             incoming_bend_end_label[side_ccw], num_incoming_bend[side_ccw],
                                                                                                             wire_bend_label[to_side], num_bend_wires[to_side]);
                            }
                            else if (seg_details[itrack].switch_dir_type() == DOWN_TYPE && is_corner(grid, i, j)) {
                                sblock_pattern[i][j][side_ccw][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                             incoming_bend_end_label[side_ccw], num_incoming_bend[side_ccw],
                                                                                                             wire_bend_label[to_side], num_bend_wires[to_side]);
                            }
                        }
                }
            }
        }

        if (incoming_wire_label[side_opp]) {
            /* Figure out the channel and segment for a certain direction */
            vert = ((side_opp == TOP) || (side_opp == BOTTOM));
            pos_dir = ((side_opp == TOP) || (side_opp == RIGHT));
            chan = (vert ? i : j);
            sb_seg = (vert ? j : i);
            seg = (pos_dir ? (sb_seg + 1) : sb_seg);

            const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg].data() : chan_details_x[seg][chan].data());
            if (seg_details[0].length() <= 0)
                continue;

            for (int itrack = 0; itrack < nodes_per_chan->max; itrack++) {
                /* not ending wire nor passing wire with sblock */
                if (incoming_wire_label[side_opp][itrack] != UN_SET) {
                    /* corner sblocks for sure have no opposite channel segments so don't care about them */
                    if (incoming_wire_label[side_opp][itrack] < num_ending_wires[side_opp]) {
                        /* If it is bend segment, use wilton way */
                        if (seg_details[itrack].bend_len() > 1) {
                            // fringe pattern
                            if (is_style_perimeter(grid, i, j, to_side, seg_details[itrack].switch_dir_type())
                                && !is_corner(grid, i, j)) {
                                // middle part of the bend segment
                                if (seg_details[itrack].part_idx() != seg_details[itrack].bend_len() - 1) {
                                    sblock_pattern[i][j][side_opp][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                                 incoming_bend_end_label[side_opp], num_incoming_bend[side_opp],
                                                                                                                 wire_bend_label[to_side], num_bend_wires[to_side]);
                                    continue;
                                }
                            }

                           
                        } 
                    }
                }
            }
        }
    }

    for (e_side side : {TOP, RIGHT, BOTTOM, LEFT}) {
        if (incoming_wire_label[side]) {
            vtr::free(incoming_wire_label[side]);
        }
        if (wire_mux_on_track[side]) {
            vtr::free(wire_mux_on_track[side]);
        }
        if (wire_bend_label[side]) {
            vtr::free(wire_bend_label[side]);
        }
        if (incoming_bend_end_label[side]) {
            vtr::free(incoming_bend_end_label[side]);
        }
    }
}

void load_sblock_pattern_lookup(const int i,
                                const int j,
                                const DeviceGrid& grid,
                                const t_chan_width* nodes_per_chan,
                                const t_chan_details& chan_details_x,
                                const t_chan_details& chan_details_y,
                                const int /*Fs*/,
                                const enum e_switch_block_type switch_block_type,
                                t_sblock_pattern& sblock_pattern) {
    /* This routine loads a lookup table for sblock topology. The lookup table is huge
     * because the sblock varies from location to location. The i, j means the owning
     * location of the sblock under investigation. */

    /* SB's have coords from (0, 0) to (grid.width()-2, grid.height()-2) */
    VTR_ASSERT(i >= 0);
    VTR_ASSERT(i <= int(grid.width()) - 2);
    VTR_ASSERT(j >= 0);
    VTR_ASSERT(j <= int(grid.height()) - 2);

    /* May 12 - 15, 2007
     *
     * I identify three types of sblocks in the chip: 1) The core sblock, whose special
     * property is that the number of muxes (and ending wires) on each side is the same (very useful
     * property, since it leads to a N-to-N assignment problem with ending wires). 2) The corner sblock
     * which is same as a L=1 core sblock with 2 sides only (again N-to-N assignment problem). 3) The
     * fringe / chip edge sblock which is most troublesome, as balance in each side of muxes is
     * attainable but balance in the entire sblock is not. The following code first identifies the
     * incoming wires, which can be classified into incoming passing wires with sblock and incoming
     * ending wires (the word "incoming" is sometimes dropped for ease of discussion). It appropriately
     * labels all the wires on each side by the following order: By the call to label_incoming_wires,
     * which labels for one side, the order is such that the incoming ending wires (always with sblock)
     * are labelled first 0,1,2,... p-1, then the incoming passing wires with sblock are labelled
     * p,p+1,p+2,... k-1 (for total of k). By this convention, one can easily distinguish the ending
     * wires from the passing wires by checking a label against num_ending_wires variable.
     *
     * After labelling all the incoming wires, this routine labels the muxes on the side we're currently
     * connecting to (iterated for four sides of the sblock), called the to_side. The label scheme is
     * the natural order of the muxes by their track #. Also we find the number of muxes.
     *
     * For each to_side, the total incoming wires that connect to the muxes on to_side
     * come from three sides: side_1 (to_side's right), side_2 (to_side's left) and opp_side.
     * The problem of balancing mux size is then: considering all incoming passing wires
     * with sblock on side_1, side_2 and opp_side, how to assign them to the muxes on to_side
     * (with specified Fs) in a way that mux size is imbalanced by at most 1. I solve this
     * problem by this approach: the first incoming passing wire will connect to 0, 1, 2,
     * ..., Fs_per_side - 1, then the next incoming passing wire will connect to
     * Fs_per_side, Fs_per_side+1, ..., Fs_per_side*2-1, and so on. This consistent STAGGERING
     * ensures N-to-N assignment is perfectly balanced and M-to-N assignment is imbalanced by no
     * more than 1.
     */

    /* SB's range from (0, 0) to (grid.width() - 2, grid.height() - 2) */
    /* First find all four sides' incoming wires */
    int* wire_mux_on_track[4];
    int* incoming_wire_label[4];
    int num_incoming_wires[4];
    int num_ending_wires[4];
    int num_wire_muxes[4];

    int* wire_bend_label[4];
    int num_bend_wires[4];
    int* incoming_bend_end_label[4];
    int num_incoming_bend[4];

    bool is_dangling_side[4];

    /* "Label" if the side of switch block is the dangling side according to the location of switch block */
    for (e_side side : {TOP, RIGHT, BOTTOM, LEFT}) {
        is_dangling_side[side] = label_dangling_side(side, i, j, grid);
    }

    /* "Label" the wires around the switch block by connectivity. */
    for (e_side side : {TOP, RIGHT, BOTTOM, LEFT}) {
        /* Assume the channel segment doesn't exist. */
        wire_mux_on_track[side] = nullptr;
        incoming_wire_label[side] = nullptr;
        num_incoming_wires[side] = 0;
        num_ending_wires[side] = 0;
        num_wire_muxes[side] = 0;

        wire_bend_label[side] = nullptr;
        num_bend_wires[side] = 0;
        incoming_bend_end_label[side] = nullptr;
        num_incoming_bend[side] = 0;

        /* Skip the side and leave the zero'd value if the
         * channel segment doesn't exist. */
        bool skip = true;
        switch (side) {
            case TOP:
                if (j < int(grid.height()) - 2) {
                    skip = false;
                }
                break;
            case RIGHT:
                if (i < int(grid.width()) - 2) {
                    skip = false;
                }
                break;
            case BOTTOM:
                if (j > 0) {
                    skip = false;
                }
                break;
            case LEFT:
                if (i > 0) {
                    skip = false;
                }
                break;
            default:
                VTR_ASSERT_MSG(false, "Unrecognzied side");
                break;
        }
        if (skip) {
            continue;
        }

        /* Figure out the channel and segment for a certain direction */
        bool vert = ((side == TOP) || (side == BOTTOM));
        bool pos_dir = ((side == TOP) || (side == RIGHT));
        int chan_len = (vert ? grid.height() : grid.width()) - 2; //-2 for no perim channels
        int chan = (vert ? i : j);
        int sb_seg = (vert ? j : i);
        int seg = (pos_dir ? (sb_seg + 1) : sb_seg);

        const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg] : chan_details_x[seg][chan]).data();
        if (seg_details[0].length() <= 0)
            continue;

        /* Figure out all the tracks on a side that are ending and the
         * ones that are passing through and have a SB. */
        enum e_direction end_dir = (pos_dir ? DEC_DIRECTION : INC_DIRECTION);
        incoming_wire_label[side] = label_incoming_wires(chan, seg, sb_seg,
                                                         seg_details, chan_len, end_dir, nodes_per_chan->max,
                                                         &num_incoming_wires[side], &num_ending_wires[side]);

        incoming_bend_end_label[side] = label_incoming_bend_wires_ending(chan, seg,
                                                                         seg_details, UNDEFINED, chan_len, end_dir, nodes_per_chan->max,
                                                                         &num_incoming_bend[side], is_dangling_side[side]);

        /* Figure out all the tracks on a side that are starting. */
        int dummy;
        enum e_direction start_dir = (pos_dir ? INC_DIRECTION : DEC_DIRECTION);
        wire_mux_on_track[side] = label_wire_muxes(chan, seg,
                                                   seg_details, UNDEFINED, chan_len, start_dir, nodes_per_chan->max,
                                                   false, &num_wire_muxes[side], &dummy, is_dangling_side[side]);

        wire_bend_label[side] = label_wire_bends(chan, seg,
                                                 seg_details, UNDEFINED, chan_len, start_dir, nodes_per_chan->max,
                                                 &num_bend_wires[side], is_dangling_side[side]);
    }

    for (e_side to_side : {TOP, RIGHT, BOTTOM, LEFT}) {
        /* Can't do anything if no muxes on this side. */
        if (num_wire_muxes[to_side] == 0)
            continue;

        /* Figure out side rotations */
        VTR_ASSERT((TOP == 0) && (RIGHT == 1) && (BOTTOM == 2) && (LEFT == 3));
        int side_cw = (to_side + 1) % 4;
        int side_opp = (to_side + 2) % 4;
        int side_ccw = (to_side + 3) % 4;

        /* Figure out the channel and segment for a certain direction */
        bool vert = ((to_side == TOP) || (to_side == BOTTOM));
        bool pos_dir = ((to_side == TOP) || (to_side == RIGHT));
        int chan = (vert ? i : j);
        int sb_seg = (vert ? j : i);
        int seg = (pos_dir ? (sb_seg + 1) : sb_seg);
        int chan_len = (vert ? grid.height() : grid.width()) - 2; //-2 for no perim channels

        const t_chan_seg_details* to_seg_details = (vert ? chan_details_y[chan][seg].data() : chan_details_x[seg][chan].data());
        if (to_seg_details[0].length() <= 0)
            continue;

        /* For the core sblock:
         * The new order for passing wires should appear as
         * 0,1,2..,scw-1, for passing wires with sblock on side_cw
         * scw,scw+1,...,sccw-1, for passing wires with sblock on side_ccw
         * sccw,sccw+1,... for passing wires with sblock on side_opp.
         * This way, I can keep the imbalance to at most 1.
         *
         * For the fringe sblocks, I don't distinguish between
         * passing and ending wires so the above statement still holds
         * if you replace "passing" by "incoming" */

        if (incoming_wire_label[side_cw]) {
            /* Figure out the channel and segment for a certain direction */
            vert = ((side_cw == TOP) || (side_cw == BOTTOM));
            pos_dir = ((side_cw == TOP) || (side_cw == RIGHT));
            chan_len = (vert ? grid.height() : grid.width()) - 2; //-2 for no perim channels
            chan = (vert ? i : j);
            sb_seg = (vert ? j : i);
            seg = (pos_dir ? (sb_seg + 1) : sb_seg);

            const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg].data() : chan_details_x[seg][chan].data());
            if (seg_details[0].length() <= 0)
                continue;

            for (int ichan = 0; ichan < nodes_per_chan->max; ichan++) { //nodes_per_chan->max == max_chan_width 不含中间节点
                int itrack = ichan;
                if (side_cw == TOP || side_cw == BOTTOM) {
                    itrack = ichan % nodes_per_chan->y_list[i];
                } else if (side_cw == RIGHT || side_cw == LEFT) {
                    itrack = ichan % nodes_per_chan->x_list[j];
                }

                /* DOWN stype from side_cw , if incoming wire is the ending of bend segment,
                 * and the bend style is the DOWN_STYLE
                 * then finding the track map bewtween side_cw and to_side
                 */

                // bend ending point
                e_direction check_dir = (pos_dir ? DEC_DIRECTION : INC_DIRECTION);
                if (incoming_wire_label[side_cw][itrack] != UN_SET
                    && incoming_wire_label[side_cw][itrack] < num_ending_wires[side_cw]
                    && seg_details[itrack].direction() == check_dir) {
                    // not the dangling incoming segment
                    int true_seg_start = seg - (seg + seg_details[itrack].length() + chan - seg_details[itrack].start()) % seg_details[itrack].length();
                    int true_seg_end = true_seg_start + seg_details[itrack].length() - 1;
                    int start = get_seg_start(seg_details, itrack, chan, seg);
                    int end = get_seg_end(seg_details, itrack, start, chan, chan_len);
                    bool clipped = (check_dir == DEC_DIRECTION ? true_seg_start != start : true_seg_end != end);
                    if ((!is_dangling_side[side_cw]) || !clipped)
                        // middle part of the bend segment
                        if (seg_details[itrack].part_idx() != seg_details[itrack].bend_len() - 1) {
                            // if bend style is the DOWN_TYPE
                            if (seg_details[itrack].switch_dir_type() == DOWN_TYPE) {
                                // assign to_track
                                sblock_pattern[i][j][side_cw][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                            incoming_bend_end_label[side_cw], num_incoming_bend[side_cw],
                                                                                                            wire_bend_label[to_side], num_bend_wires[to_side]);
                                continue;
                            }
                            if (seg_details[itrack].switch_dir_type() == UP_TYPE && is_corner(grid, i, j)) {
                                sblock_pattern[i][j][side_cw][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                            incoming_bend_end_label[side_cw], num_incoming_bend[side_cw],
                                                                                                            wire_bend_label[to_side], num_bend_wires[to_side]);
                                continue;
                            }
                        }
                }

                if (incoming_wire_label[side_cw][itrack] != UN_SET) {
                    int mux = get_simple_switch_block_track((enum e_side)side_cw,
                                                            (enum e_side)to_side,
                                                            incoming_wire_label[side_cw][ichan],
                                                            switch_block_type,
                                                            num_wire_muxes[to_side]);

                    if (sblock_pattern[i][j][side_cw][to_side][itrack][0] == UN_SET) {
                        sblock_pattern[i][j][side_cw][to_side][itrack][0] = mux;
                    } else if (sblock_pattern[i][j][side_cw][to_side][itrack][2] == UN_SET) {
                        sblock_pattern[i][j][side_cw][to_side][itrack][2] = mux;
                    }
                }
            }
        }

        if (incoming_wire_label[side_ccw]) {
            /* Figure out the channel and segment for a certain direction */
            vert = ((side_ccw == TOP) || (side_ccw == BOTTOM));
            pos_dir = ((side_ccw == TOP) || (side_ccw == RIGHT));
            chan = (vert ? i : j);
            sb_seg = (vert ? j : i);
            seg = (pos_dir ? (sb_seg + 1) : sb_seg);

            const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg].data() : chan_details_x[seg][chan].data());
            if (seg_details[0].length() <= 0)
                continue;

            for (int ichan = 0; ichan < nodes_per_chan->max; ichan++) {
                int itrack = ichan;
                if (side_ccw == TOP || side_ccw == BOTTOM) {
                    itrack = ichan % nodes_per_chan->y_list[i];
                } else if (side_ccw == RIGHT || side_ccw == LEFT) {
                    itrack = ichan % nodes_per_chan->x_list[j];
                }

                /* UP stype from side_ccw , if incoming wire is the ending of bend segment,
                 * and the bend style is the UP_STYPE
                 * then finding the track map bewtween side_ccw and to_side
                 */

                // bend ending point
                e_direction check_dir = (pos_dir ? DEC_DIRECTION : INC_DIRECTION);
                if (incoming_wire_label[side_ccw][itrack]!= UN_SET
                    && incoming_wire_label[side_ccw][itrack] < num_ending_wires[side_ccw]
                    && seg_details[itrack].direction() == check_dir) {
                    // not the dangling incoming segment
                    int true_seg_start = seg - (seg + seg_details[itrack].length() + chan - seg_details[itrack].start()) % seg_details[itrack].length();
                    int true_seg_end = true_seg_start + seg_details[itrack].length() - 1;
                    int start = get_seg_start(seg_details, itrack, chan, seg);
                    int end = get_seg_end(seg_details, itrack, start, chan, chan_len);
                    bool clipped = (check_dir == DEC_DIRECTION ? true_seg_start != start : true_seg_end != end);
                    if (!clipped || (!is_dangling_side[side_ccw]))
                        // middle part of the bend segment
                        if (seg_details[itrack].part_idx() != seg_details[itrack].bend_len() - 1) {
                            // if bend style is the UP_TYPE
                            if (seg_details[itrack].switch_dir_type() == UP_TYPE) {
                                // assign to_track
                                sblock_pattern[i][j][side_ccw][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                             incoming_bend_end_label[side_ccw], num_incoming_bend[side_ccw],
                                                                                                             wire_bend_label[to_side], num_bend_wires[to_side]);
                                continue;
                            }
                            if (seg_details[itrack].switch_dir_type() == DOWN_TYPE && is_corner(grid, i, j)) {
                                sblock_pattern[i][j][side_ccw][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                             incoming_bend_end_label[side_ccw], num_incoming_bend[side_ccw],
                                                                                                             wire_bend_label[to_side], num_bend_wires[to_side]);
                                continue;
                            }
                        }
                }

                if (incoming_wire_label[side_ccw][itrack] != UN_SET) {
                    int mux = get_simple_switch_block_track((enum e_side)side_ccw,
                                                            (enum e_side)to_side,
                                                            incoming_wire_label[side_ccw][ichan],
                                                            switch_block_type, num_wire_muxes[to_side]);

                    if (sblock_pattern[i][j][side_ccw][to_side][itrack][0] == UN_SET) {
                        sblock_pattern[i][j][side_ccw][to_side][itrack][0] = mux;
                    } else if (sblock_pattern[i][j][side_ccw][to_side][itrack][2] == UN_SET) {
                        sblock_pattern[i][j][side_ccw][to_side][itrack][2] = mux;
                    }
                }
            }
        }

        if (incoming_wire_label[side_opp]) {
            /* Figure out the channel and segment for a certain direction */
            vert = ((side_opp == TOP) || (side_opp == BOTTOM));
            pos_dir = ((side_opp == TOP) || (side_opp == RIGHT));
            chan = (vert ? i : j);
            sb_seg = (vert ? j : i);
            seg = (pos_dir ? (sb_seg + 1) : sb_seg);

            const t_chan_seg_details* seg_details = (vert ? chan_details_y[chan][seg].data() : chan_details_x[seg][chan].data());
            if (seg_details[0].length() <= 0)
                continue;

            for (int itrack = 0; itrack < nodes_per_chan->max; itrack++) {
                /* not ending wire nor passing wire with sblock */
                if (incoming_wire_label[side_opp][itrack] != UN_SET) {
                    /* corner sblocks for sure have no opposite channel segments so don't care about them */
                    if (incoming_wire_label[side_opp][itrack] < num_ending_wires[side_opp]) {
                        /* If it is bend segment, use wilton way */
                        if (seg_details[itrack].bend_len() > 1) {
                            // fringe pattern
                            if (is_style_perimeter(grid, i, j, to_side, seg_details[itrack].switch_dir_type())
                                && !is_corner(grid, i, j)) {
                                // middle part of the bend segment
                                if (seg_details[itrack].part_idx() != seg_details[itrack].bend_len() - 1) {
                                    sblock_pattern[i][j][side_opp][to_side][itrack][1] = find_to_side_bend_point(itrack, seg_details, to_seg_details,
                                                                                                                 incoming_bend_end_label[side_opp], num_incoming_bend[side_opp],
                                                                                                                 wire_bend_label[to_side], num_bend_wires[to_side]);
                                    continue;
                                }
                            }

                            //TODO: change this pattern, investigate different patterns
                            int mux = get_wilton_like_straight_sblock((enum e_side)to_side,
                                                                      incoming_wire_label[side_opp][itrack],
                                                                      seg_details[itrack].switch_dir_type(),
                                                                      num_wire_muxes[to_side]);

                            VTR_ASSERT(sblock_pattern[i][j][side_opp][to_side][itrack][0] == UN_SET);
                            sblock_pattern[i][j][side_opp][to_side][itrack][0] = mux;
                        } else {
                            /* The ending wires in core sblocks form N-to-N assignment problem, so can
                             * use any pattern such as Wilton */
                            /* In the direct connect case, I know for sure the init mux is at the same track #
                             * as this ending wire, but still need to find the init mux label for Fs > 3 */
                            int mux = find_label_of_track(wire_mux_on_track[to_side],
                                                          num_wire_muxes[to_side], itrack);
                            sblock_pattern[i][j][side_opp][to_side][itrack][0] = mux;
                        }

                    } else {
                        /* These are wire segments that pass through the switch block.
                         *
                         * There is no connection from wire segment midpoints to the opposite switch block
                         * side, so there's nothing to be done here (at Fs=3, this connection is implicit for passing
                         * wires and at Fs>3 the code in this function seems to create heavily unbalanced
                         * switch patterns). Additionally, the code in build_rr_chan() explicitly skips
                         * connections from wire segment midpoints to the opposide sb side (for switch block patterns
                         * generated with this function) so any such assignment to sblock_pattern will be ignored anyway. */
                    }
                }
            }
        }
    }

    for (e_side side : {TOP, RIGHT, BOTTOM, LEFT}) {
        if (incoming_wire_label[side]) {
            vtr::free(incoming_wire_label[side]);
        }
        if (wire_mux_on_track[side]) {
            vtr::free(wire_mux_on_track[side]);
        }
        if (wire_bend_label[side]) {
            vtr::free(wire_bend_label[side]);
        }
        if (incoming_bend_end_label[side]) {
            vtr::free(incoming_bend_end_label[side]);
        }
    }
}

/* checks if the specified coordinates represent a corner of the FPGA */
static bool is_corner(const DeviceGrid& grid, int x, int y) {
    bool is_corner = false;
    if ((x == 0 && y == 0) || (x == 0 && y == int(grid.height()) - 2) || //-2 for no perim channels
        (x == int(grid.width()) - 2 && y == 0) ||                        //-2 for no perim channels
        (x == int(grid.width()) - 2 && y == int(grid.height()) - 2)) {   //-2 for no perim channels
        is_corner = true;
    }
    return is_corner;
}

static bool is_perimeter_left(const DeviceGrid& grid, int x, int y) {
    bool is_perimeter_left = false;
    if (x == 0) {
        is_perimeter_left = true;
    }
    return is_perimeter_left;
}

static bool is_perimeter_right(const DeviceGrid& grid, int x, int y) {
    bool is_perimeter_right = false;
    if (x == int(grid.width()) - 2) {
        is_perimeter_right = true;
    }
    return is_perimeter_right;
}

static bool is_perimeter_top(const DeviceGrid& grid, int x, int y) {
    bool is_perimeter_top = false;
    if (y == int(grid.height()) - 2) {
        is_perimeter_top = true;
    }
    return is_perimeter_top;
}

static bool is_perimeter_bottom(const DeviceGrid& grid, int x, int y) {
    bool is_perimeter_bottom = false;
    if (y == 0) {
        is_perimeter_bottom = true;
    }
    return is_perimeter_bottom;
}

static bool is_style_perimeter(const DeviceGrid& grid, int x, int y, e_side to_side, e_switch_dir_type switch_type) {
    bool is_style_perimeter = false;
    if (switch_type == DOWN_TYPE) {
        if (is_perimeter_right(grid, x, y) && to_side == TOP)
            is_style_perimeter = true;
        else if (is_perimeter_left(grid, x, y) && to_side == BOTTOM)
            is_style_perimeter = true;
        else if (is_perimeter_top(grid, x, y) && to_side == LEFT)
            is_style_perimeter = true;
        else if (is_perimeter_bottom(grid, x, y) && to_side == RIGHT)
            is_style_perimeter = true;
    } else if (switch_type == UP_TYPE) {
        if (is_perimeter_right(grid, x, y) && to_side == BOTTOM)
            is_style_perimeter = true;
        else if (is_perimeter_left(grid, x, y) && to_side == TOP)
            is_style_perimeter = true;
        else if (is_perimeter_top(grid, x, y) && to_side == RIGHT)
            is_style_perimeter = true;
        else if (is_perimeter_bottom(grid, x, y) && to_side == LEFT)
            is_style_perimeter = true;
    }
    return is_style_perimeter;
}

static bool label_dangling_side(const e_side side, const int x, const int y, const DeviceGrid& grid) {
    bool is_dangling = false;
    // left-bottom corner
    if (x == 0 && y == 0) {
        if (side == TOP || side == RIGHT)
            is_dangling = true;
    }
    // left-top corner
    else if (x == 0 && y == int(grid.height()) - 2) {
        if (side == BOTTOM || side == RIGHT)
            is_dangling = true;
    }
    // right-bottom corner
    else if (x == int(grid.width()) - 2 && y == 0) {
        if (side == TOP || side == LEFT)
            is_dangling = true;
    }
    // right-top corner
    else if (x == int(grid.width()) - 2 && y == int(grid.height()) - 2) {
        if (side == BOTTOM || side == LEFT)
            is_dangling = true;
    }
    // top fringe
    else if (y == int(grid.height()) - 2) {
        if (side == BOTTOM)
            is_dangling = true;
    }
    // bottom fringe
    else if (y == 0) {
        if (side == TOP)
            is_dangling = true;
    }
    // right fringe
    else if (x == int(grid.width()) - 2) {
        if (side == LEFT)
            is_dangling = true;
    }
    // left fringe
    else if (x == 0) {
        if (side == RIGHT)
            is_dangling = true;
    }

    return is_dangling;
}

static int* label_wire_bends(const int chan_num,
                             const int seg_num,
                             const t_chan_seg_details* seg_details,
                             const int seg_type_index,
                             const int max_len,
                             const enum e_direction dir,
                             const int max_chan_width,
                             int* num_bend_wires,
                             bool is_dangling_side) {
    /* Labels the bend segment start points on that side (seg_num, chan_num, direction). The returned array
     * maps a label to the actual track #: array[0] = <the track number of the first/lowest mux>
     * This routine orders wire bends by their natural order, i.e. track #
     * If seg_type_index == UNDEFINED, all segments in the channel are considered. Otherwise this routine
     * only looks at segments that belong to the specified segment type. */

    int itrack, start, end, pass;
    int* labels = nullptr;
    bool is_endpoint, isdangling;
    int num_labels;
    int true_seg_start, true_seg_end, length;
    bool clipped;

    /* COUNT pass then a LOAD pass */
    num_labels = 0;
    for (pass = 0; pass < 2; ++pass) {
        /* Alloc the list on LOAD pass */
        if (pass > 0) {
            labels = (int*)vtr::malloc(sizeof(int) * num_labels);
            num_labels = 0;
        }

        /* Find the tracks that are starting. */
        for (itrack = 0; itrack < max_chan_width; ++itrack) {
            start = get_seg_start(seg_details, itrack, chan_num, seg_num);
            end = get_seg_end(seg_details, itrack, start, chan_num, max_len);

            /* Skip tracks that are undefined */
            if (seg_details[itrack].length() == 0) {
                continue;
            }

            /* Skip tracks going the wrong way */
            if (seg_details[itrack].direction() != dir) {
                continue;
            }

            if (seg_type_index != UNDEFINED) {
                /* skip tracks that don't belong to the specified segment type */
                if (seg_details[itrack].index() != seg_type_index) {
                    continue;
                }
            }

            length = seg_details[itrack].length();
            true_seg_start = seg_num - (seg_num + length + chan_num - seg_details[itrack].start()) % length;
            true_seg_end = true_seg_start + seg_details[itrack].length() - 1;
            clipped = (dir == INC_DIRECTION ? true_seg_start != start : true_seg_end != end);

            /* Determine if we are a wire startpoint */
            is_endpoint = (seg_num == start);
            if (DEC_DIRECTION == seg_details[itrack].direction()) {
                is_endpoint = (seg_num == end);
            }

            isdangling = false;
            if (is_endpoint && is_dangling_side) {
                if (clipped) {
                    isdangling = true;
                }
            }

            /* Skip tracks that in the middle of bend segment */
            if (seg_details[itrack].part_idx() > 0 && is_endpoint && !isdangling) {
                if (pass > 0) {
                    labels[num_labels] = itrack;
                }
                num_labels++;
                continue;
            }
        }
    }

    *num_bend_wires = num_labels;

    return labels;
}

static int* label_wire_muxes(const int chan_num,
                             const int seg_num,
                             const t_chan_seg_details* seg_details,
                             const int seg_type_index,
                             const int max_len,
                             const enum e_direction dir,
                             const int max_chan_width,
                             const bool check_cb,
                             int* num_wire_muxes,
                             int* num_wire_muxes_cb_restricted,
                             bool is_dangling_side) {
    /* Labels the muxes on that side (seg_num, chan_num, direction). The returned array
     * maps a label to the actual track #: array[0] = <the track number of the first/lowest mux>
     * This routine orders wire muxes by their natural order, i.e. track #
     * If seg_type_index == UNDEFINED, all segments in the channel are considered. Otherwise this routine
     * only looks at segments that belong to the specified segment type. */

    int itrack, start, end, num_labels, num_labels_restricted, pass;
    int* labels = nullptr;
    bool is_endpoint, isdangling;
    int true_seg_start, true_seg_end, length;
    bool clipped;

    /* COUNT pass then a LOAD pass */
    num_labels = 0;
    num_labels_restricted = 0;
    for (pass = 0; pass < 2; ++pass) {
        /* Alloc the list on LOAD pass */
        if (pass > 0) {
            labels = (int*)vtr::malloc(sizeof(int) * num_labels);
            num_labels = 0;
        }

        /* Find the tracks that are starting. */
        for (itrack = 0; itrack < max_chan_width; ++itrack) {
            start = get_seg_start(seg_details, itrack, chan_num, seg_num);
            end = get_seg_end(seg_details, itrack, start, chan_num, max_len);

            /* Skip tracks that are undefined */
            if (seg_details[itrack].length() == 0) {
                continue;
            }

            /* Skip tracks going the wrong way */
            if (seg_details[itrack].direction() != dir) {
                continue;
            }

            if (seg_type_index != UNDEFINED) {
                /* skip tracks that don't belong to the specified segment type */
                if (seg_details[itrack].index() != seg_type_index) {
                    continue;
                }
            }

            length = seg_details[itrack].length();
            true_seg_start = seg_num - (seg_num + length + chan_num - seg_details[itrack].start()) % length;
            true_seg_end = true_seg_start + seg_details[itrack].length() - 1;
            clipped = (dir == INC_DIRECTION ? true_seg_start != start : true_seg_end != end);

            /* Determine if we are a wire startpoint */
            is_endpoint = (seg_num == start);
            if (DEC_DIRECTION == seg_details[itrack].direction()) {
                is_endpoint = (seg_num == end);
            }

            isdangling = false;
            if (is_endpoint && is_dangling_side) {
                if (clipped) {
                    isdangling = true;
                }
            }

            if (seg_details[itrack].part_idx() > 0 && !isdangling)
                continue;
            /* Count the labels and load if LOAD pass */
            if (is_endpoint) {
                /*
                 * not all wire endpoints can be driven by OPIN (depending on the <cb> pattern in the arch file)
                 * the check_cb is targeting this arch specification:
                 * if this function is called by get_unidir_opin_connections(),
                 * then we need to check if mux connections can be added to this type of wire,
                 * otherwise, this function should not consider <cb> specification.
                 */
                if ((!check_cb) || (seg_details[itrack].cb(0) == true)) {
                    if (pass > 0) {
                        labels[num_labels] = itrack;
                    }
                    ++num_labels;
                }
                if (pass > 0)
                    num_labels_restricted += (seg_details[itrack].cb(0) == true) ? 1 : 0;
            }
        }
    }

    *num_wire_muxes = num_labels;
    *num_wire_muxes_cb_restricted = num_labels_restricted;

    return labels;
}

static int* label_incoming_bend_wires_ending(const int chan_num,
                                             const int seg_num,
                                             const t_chan_seg_details* seg_details,
                                             const int seg_type_index,
                                             const int max_len,
                                             const enum e_direction dir,
                                             const int max_chan_width,
                                             int* num_incoming_bend,
                                             bool is_dangling_side) {
    /* Labels the incoming ending of bend wires on that side (seg_num, chan_num, direction).
     * The returned array maps a track # to a label: array[0] = <the new hash value/label for track 0>,
     * the labels 0,1,2,.. identify consecutive incoming wires that have sblock (passing wires with sblock and ending wires) */

    int itrack, start, end, pass;
    int* labels = nullptr;
    bool is_endpoint, isdangling;
    int num_labels;
    int true_seg_start, true_seg_end, length;
    bool clipped;

    /* COUNT pass then a LOAD pass */
    num_labels = 0;
    for (pass = 0; pass < 2; ++pass) {
        /* Alloc the list on LOAD pass */
        if (pass > 0) {
            labels = (int*)vtr::malloc(sizeof(int) * num_labels);
            num_labels = 0;
        }

        /* Find the tracks that are starting. */
        for (itrack = 0; itrack < max_chan_width; ++itrack) {
            start = get_seg_start(seg_details, itrack, chan_num, seg_num);
            end = get_seg_end(seg_details, itrack, start, chan_num, max_len);

            /* Skip tracks that are undefined */
            if (seg_details[itrack].length() == 0) {
                continue;
            }

            /* Skip tracks going the wrong way */
            if (seg_details[itrack].direction() != dir) {
                continue;
            }

            if (seg_type_index != UNDEFINED) {
                /* skip tracks that don't belong to the specified segment type */
                if (seg_details[itrack].index() != seg_type_index) {
                    continue;
                }
            }

            length = seg_details[itrack].length();
            true_seg_start = seg_num - (seg_num + length + chan_num - seg_details[itrack].start()) % length;
            true_seg_end = true_seg_start + seg_details[itrack].length() - 1;
            clipped = (dir == DEC_DIRECTION ? true_seg_start != start : true_seg_end != end);

            /* Determine if we are a wire endpoint */
            is_endpoint = (seg_num == end);
            if (DEC_DIRECTION == seg_details[itrack].direction()) {
                is_endpoint = (seg_num == start);
            }

            isdangling = false;
            if (is_endpoint && is_dangling_side) {
                if (clipped) {
                    isdangling = true;
                }
            }

            /* Skip tracks that in the middle of bend segment */
            if (seg_details[itrack].part_idx() < seg_details[itrack].bend_len() - 1 && is_endpoint && !isdangling) {
                if (pass > 0) {
                    labels[num_labels] = itrack;
                }
                num_labels++;
                continue;
            }
        }
    }

    *num_incoming_bend = num_labels;

    return labels;
}

static int* label_incoming_wires(const int chan_num,
                                 const int seg_num,
                                 const int sb_seg,
                                 const t_chan_seg_details* seg_details,
                                 const int max_len,
                                 const enum e_direction dir,
                                 const int max_chan_width,
                                 int* num_incoming_wires,
                                 int* num_ending_wires) {
    /* Labels the incoming wires on that side (seg_num, chan_num, direction).
     * The returned array maps a track # to a label: array[0] = <the new hash value/label for track 0>,
     * the labels 0,1,2,.. identify consecutive incoming wires that have sblock (passing wires with sblock and ending wires) */

    int itrack, start, end, i, num_passing, num_ending, pass;
    int* labels;
    bool sblock_exists, is_endpoint;

    /* Alloc the list of labels for the tracks */
    labels = (int*)vtr::malloc(max_chan_width * sizeof(int));
    for (i = 0; i < max_chan_width; ++i) {
        labels[i] = UN_SET; /* crash hard if unset */
    }

    num_ending = 0;
    num_passing = 0;
    for (pass = 0; pass < 2; ++pass) {
        for (itrack = 0; itrack < max_chan_width; ++itrack) {
            /* Skip tracks that are undefined */
            if (seg_details[itrack].length() == 0) {
                continue;
            }

            if (seg_details[itrack].direction() == dir) {
                start = get_seg_start(seg_details, itrack, chan_num, seg_num);
                end = get_seg_end(seg_details, itrack, start, chan_num, max_len);

                /* Determine if we are a wire endpoint */
                is_endpoint = (seg_num == end);
                if (DEC_DIRECTION == seg_details[itrack].direction()) {
                    is_endpoint = (seg_num == start);
                }

                /* Determine if we have a sblock on the wire */
                sblock_exists = is_sblock(chan_num, seg_num, sb_seg, itrack,
                                          seg_details, UNI_DIRECTIONAL);

                switch (pass) {
                        /* On first pass, only load ending wire labels. */
                    case 0:
                        if (is_endpoint) {
                            labels[itrack] = num_ending;
                            ++num_ending;
                        }
                        break;

                        /* On second pass, load the passing wire labels. They
                         * will follow after the ending wire labels. */
                    case 1:
                        if ((false == is_endpoint) && sblock_exists) {
                            labels[itrack] = num_ending + num_passing;
                            ++num_passing;
                        }
                        break;
                    default:
                        VTR_ASSERT_MSG(false, "Unrecognized pass");
                        break;
                }
            }
        }
    }

    *num_incoming_wires = num_passing + num_ending;
    *num_ending_wires = num_ending;
    return labels;
}

static int find_label_of_track(int* wire_mux_on_track,
                               int num_wire_muxes,
                               int from_track) {
    /* Returns the index/label in array wire_mux_on_track whose entry equals from_track. If none are
     * found, then returns the index of the entry whose value is the largest */
    int i_label = -1;
    int max_track = -1;

    for (int i = 0; i < num_wire_muxes; i++) {
        if (wire_mux_on_track[i] == from_track) {
            i_label = i;
            break;
        } else if (wire_mux_on_track[i] > max_track) {
            i_label = i;
            max_track = wire_mux_on_track[i];
        }
    }
    return i_label;
}

static int should_create_switchblock(const DeviceGrid& grid, int from_chan_coord, int from_seg_coord, t_rr_type from_chan_type, t_rr_type to_chan_type) {
    //Convert the chan/seg indicies to real x/y coordinates
    int y_coord;
    int x_coord;
    if (from_chan_type == CHANX) {
        y_coord = from_chan_coord;
        x_coord = from_seg_coord;
    } else {
        VTR_ASSERT(from_chan_type == CHANY);
        y_coord = from_seg_coord;
        x_coord = from_chan_coord;
    }

    t_type_ptr blk_type = grid[x_coord][y_coord].type;
    int width_offset = grid[x_coord][y_coord].width_offset;
    int height_offset = grid[x_coord][y_coord].height_offset;

    e_sb_type sb_type = blk_type->switchblock_locations[width_offset][height_offset];
    auto switch_override = blk_type->switchblock_switch_overrides[width_offset][height_offset];

    if (sb_type == e_sb_type::FULL) {
        return switch_override;
    } else if (sb_type == e_sb_type::STRAIGHT && from_chan_type == to_chan_type) {
        return switch_override;
    } else if (sb_type == e_sb_type::TURNS && from_chan_type != to_chan_type) {
        return switch_override;
    } else if (sb_type == e_sb_type::HORIZONTAL && from_chan_type == CHANX && to_chan_type == CHANX) {
        return switch_override;
    } else if (sb_type == e_sb_type::VERTICAL && from_chan_type == CHANY && to_chan_type == CHANY) {
        return switch_override;
    }

    return NO_SWITCH;
}

static bool should_apply_switch_override(int switch_override) {
    if (switch_override != NO_SWITCH && switch_override != DEFAULT_SWITCH) {
        VTR_ASSERT(switch_override >= 0);
        return true;
    }
    return false;
}

void partition_rr_graph_edges(DeviceContext& device_ctx) {
    for (size_t inode = 0; inode < device_ctx.rr_nodes.size(); ++inode) {
        device_ctx.rr_nodes[inode].partition_edges();

        VTR_ASSERT_SAFE(device_ctx.rr_nodes[inode].validate());
    }
}







int count_gsb_medium_num_v2(const t_gsb_inf& gsb, std::map<int, int>& clbpin_medium_map, std::map<std::string, int>& mux_name_medium_map) {
    auto& device_ctx = g_vpr_ctx.device();
    int i, j, k;
    int totals = 0;
    clbpin_medium_map.clear();

    for (i = 0; i < device_ctx.num_block_types; i++) {
        if (strcmp(device_ctx.block_types[i].name, gsb.pbtype_names[0].c_str()) == 0) {
            break;
        }
    }
    if (i >= device_ctx.num_block_types) {
        vpr_throw(VPR_ERROR_ARCH, get_arch_file_name(), 1, "Unable to find block %s for imux.\n", gsb.pbtype_names[0].c_str());
    }
    auto const& clb_type = device_ctx.block_types[i];
    auto const pb_type = clb_type.pb_type;

    for (i = 0; i < pb_type->num_ports; i++){
        if(pb_type->ports[i].type == IN_PORT){
            for (k = 0; k < pb_type->ports[i].num_pins; k++) {
                int pin_index;
                get_blk_pin_from_port_pin(clb_type.index, i, k, &pin_index);
                clbpin_medium_map[pin_index] = k + totals;
            }
            totals += pb_type->ports[i].num_pins;
        }
    }

    int first_stage_num = gsb.first_stages.size();

    for (i = 0; i < first_stage_num; i++) {
        mux_name_medium_map[gsb.first_stages[i].mux_name] = i + totals;
        
    }
    totals += first_stage_num;
    
    return totals;
}



void load_gsb_medium_rr_indices(const DeviceGrid& grid,
                                  const t_arch* arch,
                                  const t_medium_seg_inf& medium_seg_inf,
                                  t_medium_rr_node_indices& indices,
                                  const t_chan_details& chan_details_x,
                                  const t_chan_details& chan_details_y,
                                  int* index){
    int medium_num_xory = medium_seg_inf.gsb_mdeium_num_xory;
    auto const& L_rr_node_indices = g_vpr_ctx.mutable_device().rr_node_indices;
    for(int gsb_types_index = 0; gsb_types_index < arch->gsb_inf.size(); gsb_types_index++){
        auto const& gsb = arch->gsb_inf[gsb_types_index];
        VTR_ASSERT(indices.size() == size_t(grid.width()));
        int gsb_mdeium_num = medium_seg_inf.gsb_medium_num;
        for (int x = 0; x < grid.width()-1; x++) {
            VTR_ASSERT(indices[x].size() == size_t(grid.height()));
            for (int y = 0; y < grid.height()-1; y++) {
                VTR_ASSERT(indices[x][y].size() == size_t(1));
                VTR_ASSERT(indices[x][y][0].size() == size_t(arch->gsb_inf.size()));
                //if (grid[x][y].width_offset == 0 && grid[x][y].height_offset == 0) {//大块如何处理
                    t_type_ptr type = grid[x][y].type;
                    bool is_gsb_arch = false;
                    for (auto pb_name : gsb.pbtype_names) {
                        if (0 == strcmp(type->name, pb_name.c_str())) {
                            is_gsb_arch = true;
                            break;
                        }
                    }

                    if (is_gsb_arch) {
                        int i_gsb = 0;
                        indices[x][y][0][gsb_types_index].resize(gsb_mdeium_num, OPEN);
                        const t_chan_seg_details* seg_details_x = chan_details_x[x][y].data();
                        for (unsigned track = L_rr_node_indices[CHANX][y][x][0].size(); track < L_rr_node_indices[CHANX][y][x][0].size() + medium_num_xory; ++track) {
                            if (seg_details_x[track].is_gsb_medium()) {
                                VTR_ASSERT(seg_details_x[track].length() == 1);
                                if (i_gsb < gsb_mdeium_num) {
                                    int inode = *index;
                                    ++(*index);
                                    indices[x][y][0][gsb_types_index][i_gsb] = inode; 
                                    i_gsb++;
                                } else {
                                    break;
                                }
                            }
                        }

                        const t_chan_seg_details* seg_details_y = chan_details_y[x][y].data();
                        for (unsigned track = L_rr_node_indices[CHANY][x][y][0].size(); track < L_rr_node_indices[CHANY][x][y][0].size() + medium_num_xory; ++track) {
                            if (seg_details_y[track].is_gsb_medium()) {
                                VTR_ASSERT(seg_details_x[track].length() == 1);
                                if (i_gsb < gsb_mdeium_num) {
                                    int inode = *index;
                                    ++(*index);
                                    indices[x][y][0][gsb_types_index][i_gsb] = inode;
                                    i_gsb++;
                                } else {
                                    break;
                                }
                            }
                        }
                        //}
                    }
            }
        }
    }
}

t_medium_rr_node_indices alloc_and_load_medium_rr_node_indices(const DeviceGrid& grid,
                                                                int* index,
                                                                const t_arch* arch,
                                                                const t_medium_seg_inf& medium_seg_inf,
                                                                const t_chan_details& chan_details_x,
                                                                const t_chan_details& chan_details_y) {

    t_medium_rr_node_indices indices;//[0..grid_width-1][0..grid_height-1][OMUX/IMUX/GSB][0..num_omuxs/imuxs_types/gsb_types-1][0..pin_num-1]

    indices.resize(grid.width());
    for (size_t x = 0; x < grid.width(); ++x) {
        indices[x].resize(grid.height());
        for (size_t y = 0; y < grid.height(); ++y) {
            indices[x][y].resize(1);//imux\omux\gsb(0/1)
            indices[x][y][0].resize(arch->gsb_inf.size());
        }
    }

    load_gsb_medium_rr_indices(grid, arch, medium_seg_inf, indices, chan_details_x, chan_details_y, index);
    
    return indices;
                                                                                                                                                                                                                                                  
}

int get_track_to_track_for_gsb(int x_coord, int y_coord, const DeviceGrid& grid,
                                        const t_chan_details& chan_details_x,
                                        const t_chan_details& chan_details_y,
                                        const t_rr_node_indices& L_rr_node_indices,
                                        t_rr_edge_info_set& rr_edges_to_create,
                                        t_gsb_connection_map* gsb_conn_map,
                                        t_gsb_multistage_mux_map* gsb_multistage_mux_map){    

    e_from_type from_type = e_from_type::FT_SEGMENT;

    int edge_count = 0;
    int from_x, from_y;                     /* index into source channel */
    int to_x, to_y;                         /* index into destination channel */
    t_rr_type from_chan_type, to_chan_type; /* the type of channel - i.e. CHANX or CHANY */
    from_x = from_y = to_x = to_y = UNDEFINED;

    GSB_Lookup gsb_ms(x_coord, y_coord, e_from_type::NUM_FT_TYPES, NUM_SIDES, NUM_SIDES);
    std::unordered_map<int, std::pair<int, std::vector<int>>> to_node_loc_map;
    auto it = (*gsb_multistage_mux_map).find(gsb_ms);
    if(it != (*gsb_multistage_mux_map).end()){
        to_node_loc_map = (*it).second.to_node_loc_map;
    }

    for (e_side from_side : {TOP, RIGHT, BOTTOM, LEFT}) {
        for (e_side to_side : {TOP, RIGHT, BOTTOM, LEFT}) {
            GSB_Lookup gsb_conn(x_coord, y_coord, from_type, from_side, to_side); /* for indexing into FPGA's switchblock map */

            if (gsb_conn_map->count(gsb_conn) > 0){
                std::vector<t_gsb_edge>& conn_vector = (*gsb_conn_map)[gsb_conn];

                //to_x, to_y and so on are used, don't comment
                const t_chan_details& from_chan_details = index_into_correct_chan(x_coord, y_coord, from_side, chan_details_x, chan_details_y,
                                                                              &from_x, &from_y, &from_chan_type);
                const t_chan_details& to_chan_details = index_into_correct_chan(x_coord, y_coord, to_side, chan_details_x, chan_details_y,
                                                                        &to_x, &to_y, &to_chan_type);
                /* make sure from_x/y and to_x/y aren't out of bounds */
                if (coords_out_of_bounds(grid, to_x, to_y, to_chan_type) || coords_out_of_bounds(grid, from_x, from_y, from_chan_type)) 
                    continue;

                /* go through the connections... */
                for (int iconn = 0; iconn < (int)conn_vector.size(); ++iconn) {
                    int from_wire = conn_vector.at(iconn).from_wireOrpin;
                    int to_wire = conn_vector.at(iconn).to_wireOrpin;
                    int to_node = get_rr_node_index(L_rr_node_indices, to_x, to_y, to_chan_type, to_wire);

                    if (to_node == OPEN || to_node_loc_map.find(to_node) != to_node_loc_map.end()){
                        continue;
                    }
                    
                    int from_node = get_rr_node_index(L_rr_node_indices, from_x, from_y, from_chan_type, from_wire);

                    if (from_node == OPEN) {
                        continue;
                    }

                    int src_switch = conn_vector[iconn].switch_ind;
                    rr_edges_to_create.emplace_back(from_node, to_node, src_switch);
                    ++edge_count;
                }
            
            }
        
        }
    }
    return edge_count;
}





int get_track_to_imux_for_imux(int x_coord, int y_coord, const DeviceGrid& grid,
                                    const t_chan_details& chan_details_x,
                                    const t_chan_details& chan_details_y,
                                    const t_rr_node_indices& L_rr_node_indices,
                                    const t_medium_rr_node_indices& L_medium_rr_node_indices,
                                    t_rr_edge_info_set& rr_edges_to_create,
                                    t_gsb_connection_map* imux_conn_map,
                                    const std::map<int, int>& clbpin_to_medium,
                                    const std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map){    

    e_from_type from_type = e_from_type::FT_SEGMENT;
    int edge_count = 0;
    int from_x, from_y;                     /* index into source channel */
    t_rr_type from_chan_type; /* the type of channel - i.e. CHANX or CHANY */
    from_x = from_y = UNDEFINED;

    bool is_other_tile = false;
    auto type = grid[x_coord][y_coord].type;
    std::string pbtype_name = type->name;
    auto iter = other_pbpin_to_clbpin_map.find(pbtype_name);
    if (iter != other_pbpin_to_clbpin_map.end()) {
        is_other_tile = true;
    }

    for (e_side from_side : {TOP, RIGHT, BOTTOM, LEFT}) {
        e_side to_side = NUM_SIDES;
        GSB_Lookup gsb_conn(x_coord, y_coord, from_type, from_side, to_side); /* for indexing into FPGA's switchblock map */
        if (imux_conn_map->count(gsb_conn) > 0){
            std::vector<t_gsb_edge>& conn_vector = (*imux_conn_map)[gsb_conn];
            const t_chan_details& from_chan_details = index_into_correct_chan(x_coord, y_coord, from_side, chan_details_x, chan_details_y,
                                                                              &from_x, &from_y, &from_chan_type);
            /* make sure from_x/y and to_x/y aren't out of bounds */
            if (coords_out_of_bounds(grid, from_x, from_y, from_chan_type)) 
                continue;

            /* go through the connections... */
            for (int iconn = 0; iconn < (int)conn_vector.size(); ++iconn) {
                int from_wire = conn_vector.at(iconn).from_wireOrpin;
                int to_pb_pin = conn_vector.at(iconn).to_wireOrpin;
                int from_node = get_rr_node_index(L_rr_node_indices, from_x, from_y, from_chan_type, from_wire);
                int to_node;

                auto it = clbpin_to_medium.find(to_pb_pin);
                if(it != clbpin_to_medium.end()){
                    to_node = get_medium_rr_node_index(L_medium_rr_node_indices, x_coord, y_coord, 0, 0, it->second);//0-IMUX, 1-OMUX, 2-GSB
                }
                else {
                    if (is_other_tile) {
                        int index = grid[x_coord][y_coord].width_offset * type->height + grid[x_coord][y_coord].height_offset;
                        auto const pbpin_to_clbpin_map_tile = (iter->second)[index];
                        auto iter2 = pbpin_to_clbpin_map_tile.find(to_pb_pin);
                        if (iter2 != pbpin_to_clbpin_map_tile.end()) {
                            to_pb_pin = iter2->second;
                        } else {
                            continue;
                        }
                    }

                    for(auto side : SIDES){
                        to_node = get_rr_node_index(L_rr_node_indices, x_coord, y_coord, IPIN, to_pb_pin, side);
                        if(to_node != OPEN)
                            break;
                    }
                    //to_node = get_rr_node_index(L_rr_node_indices, x_coord, y_coord, IPIN, to_pb_pin);
                }

                if (from_node == OPEN || to_node == OPEN) {
                    continue;
                }

                int src_switch = conn_vector[iconn].switch_ind;
                rr_edges_to_create.emplace_back(from_node, to_node, src_switch);
                ++edge_count;
            }
        
        }
    }
    return edge_count;
}




int get_medium_to_pin_for_imux(int x_coord, int y_coord,
                                  const std::map<std::string, std::pair<int, int>>& record_imux_coord,
                                  const DeviceGrid& grid,
                                  const t_rr_node_indices& L_rr_node_indices,
                                  const t_medium_rr_node_indices& L_medium_rr_node_indices,
                                  t_rr_edge_info_set& rr_edges_to_create,
                                  const std::map<int, int>& clbpin_to_medium,
                                  const std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map,
                                  std::string actual_pbtype_name,
                                  const int delayless_switch){
    bool is_other_tile = false;
    auto type = grid[x_coord][y_coord].type;
    std::string pbtype_name = type->name;
    auto iter = other_pbpin_to_clbpin_map.find(pbtype_name);
    if (iter != other_pbpin_to_clbpin_map.end()) {
        is_other_tile = true;
    }

    int edge_count = 0;

    if( 0 != strcmp(pbtype_name.c_str(), actual_pbtype_name.c_str()) && !is_other_tile)
        return edge_count;

    auto it2 = record_imux_coord.find(actual_pbtype_name);
    if(it2 == record_imux_coord.end()) return edge_count;

    for(auto it : clbpin_to_medium){
        int medium_ind = it.second;
        int clb_pin = it.first;
        int from_node = get_medium_rr_node_index(L_medium_rr_node_indices, x_coord, y_coord, 0, 0, medium_ind);
        int to_node;
        if (is_other_tile) {
            int index = grid[x_coord][y_coord].width_offset * type->height + grid[x_coord][y_coord].height_offset;
            auto const pbpin_to_clbpin_map_tile = (iter->second)[index];
            auto iter2 = pbpin_to_clbpin_map_tile.find(clb_pin);
            if (iter2 != pbpin_to_clbpin_map_tile.end()) {
                clb_pin = iter2->second;
            }else{
                continue;
            }
        }
        for(auto side : SIDES){
            to_node = get_rr_node_index(L_rr_node_indices, x_coord, y_coord, IPIN, clb_pin, side);
            if(to_node != OPEN)
                break;
        }

        if (from_node == OPEN || to_node == OPEN) {
            continue;
        }

        int src_switch = delayless_switch;
        rr_edges_to_create.emplace_back(from_node, to_node, src_switch);
        ++edge_count;
    }

    return edge_count;
}




void build_other_pbpin_to_clbpin(const t_gsb_inf& gsb, std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map){
    other_pbpin_to_clbpin_map.clear();
    auto actual_pbtype = find_pb_type(gsb.pbtype_names[0]);
    VTR_ASSERT(actual_pbtype->width == 1 && actual_pbtype->height == 1);
    int actual_pbtype_ipin_num = 0;
    int actual_pbtype_opin_num = 0;
    std::vector<int> actual_pbtype_ipins;
    std::vector<int> actual_pbtype_opins;
    for (int ipin = 0; ipin < actual_pbtype->num_pins; ++ipin) {
        int iclass = actual_pbtype->pin_class[ipin];
        auto class_type = actual_pbtype->class_inf[iclass].type;

        if (class_type == DRIVER) {
            actual_pbtype_opins.push_back(ipin);
            actual_pbtype_opin_num++;
        } else {
            VTR_ASSERT(class_type == RECEIVER);
            actual_pbtype_ipins.push_back(ipin);
            actual_pbtype_ipin_num++;
        }
    }
    
    for (int i = 1; i < gsb.pbtype_names.size(); i++){
        std::vector<std::map<int, int>>& pbpin_to_clbpin = other_pbpin_to_clbpin_map[gsb.pbtype_names[i]];
        auto pb_type = find_pb_type(gsb.pbtype_names[i]);
        pbpin_to_clbpin.resize(pb_type->width * pb_type->height);
        for (int width_offset = 0; width_offset < pb_type->width; ++width_offset) {
            for (int height_offset = 0; height_offset < pb_type->height; ++height_offset) {
                std::map<int, int> pbpin_to_clbpin_tile;
                int actual_pbtype_ipin_ind = 0;
                int actual_pbtype_opin_ind = 0;
                for (int ipin = 0; ipin < pb_type->num_pins; ++ipin) {
                    for (e_side side : SIDES) {
                        if (pb_type->pinloc[width_offset][height_offset][side][ipin]) {
                            int iclass = pb_type->pin_class[ipin];
                            auto class_type = pb_type->class_inf[iclass].type;

                            if (class_type == DRIVER && actual_pbtype_opin_ind < actual_pbtype_opin_num) {
                                pbpin_to_clbpin_tile[actual_pbtype_opins[actual_pbtype_opin_ind++]] = ipin;
                            } 
                            else if(actual_pbtype_ipin_ind < actual_pbtype_ipin_num){
                                VTR_ASSERT(class_type == RECEIVER);
                                pbpin_to_clbpin_tile[actual_pbtype_ipins[actual_pbtype_ipin_ind++]] = ipin;
                            }
                        }
                    }
                }
                pbpin_to_clbpin[width_offset * pb_type->height + height_offset] = pbpin_to_clbpin_tile;
            }
        }
    }
}

int find_average_medium_rr_node_index(int device_width,
                                        int device_height,
                                        int medium_type,
                                        int ptc,
                                        const t_medium_rr_node_indices& L_medium_rr_node_indices) {
    /* Find and return the index to a rr_node that is located at the "center" *
     * of the current grid array, if possible.  In the event the "center" of  *
     * the grid array is an EMPTY or IO node, then retry alterate locations.  *
     * Worst case, this function will simply return the 1st non-EMPTY and     *
     * non-IO node.                                                           */

    int inode = get_medium_rr_node_index(L_medium_rr_node_indices, (device_width) / 2, (device_height) / 2, medium_type, 0, ptc);//默认只有一种gsb，因此num_type_index取0

    if (inode == OPEN) {
        inode = get_medium_rr_node_index(L_medium_rr_node_indices, (device_width) / 4, (device_height) / 4, medium_type, 0, ptc);
    }
    if (inode == OPEN) {
        inode = get_medium_rr_node_index(L_medium_rr_node_indices, (device_width) / 4 * 3, (device_height) / 4 * 3, medium_type, 0, ptc);
    }
    if (inode == OPEN) {
        auto& device_ctx = g_vpr_ctx.device();

        for (int x = 0; x < device_width; ++x) {
            for (int y = 0; y < device_height; ++y) {
                if (device_ctx.grid[x][y].type == device_ctx.EMPTY_TYPE)
                    continue;
                if (is_io_type(device_ctx.grid[x][y].type))
                    continue;

                inode = get_medium_rr_node_index(L_medium_rr_node_indices, x, y, medium_type, 0, ptc);
                if (inode != OPEN)
                    break;
            }
            if (inode != OPEN)
                break;
        }
    }
    return (inode);
}

/* Counts the number of wires in each wire type in the specified channel */
void count_bend_wire_type_sizes(const t_chan_seg_details* channel, int bend_start, int bend_end, std::vector<Wire_Info>& real_bend_wire_type) {
    int new_length, length;
    int new_start, start;
    int new_part, part;
    int num_wires = 0;
    Wire_Info wire_info;

    auto& device_ctx = g_vpr_ctx.device();
    auto const& segment_inf = device_ctx.arch->Segments;

    length = channel[bend_start].length();
    start = bend_start;
    part = channel[bend_start].part_idx();

    real_bend_wire_type.clear();
    for (int iwire = bend_start; iwire < bend_end; iwire++) {
        new_length = channel[iwire].length();//for bent wire we count its part length
        new_start = iwire;
        new_part = channel[iwire].part_idx();
        if (new_part != part) {
            wire_info.set(length, num_wires, start);
            real_bend_wire_type.push_back(wire_info);
            part = new_part;
            length = new_length;
            start = new_start;
            num_wires = 0;
        }
        num_wires++;
    }
    wire_info.set(length, num_wires, start);
    real_bend_wire_type.push_back(wire_info);

    return;
}

static int build_a_bound_wire(int x_coord, int y_coord,
                              const t_wire_type_sizes& wire_type_sizes,
                              const t_chan_details& chan_details,
                              const t_rr_node_indices& L_rr_node_indices,
                              t_rr_edge_info_set& rr_edges_to_create,
                              e_direction to_direction,
                              t_rr_type chan_type,
                              const int delayless_switch){
    int from_node, to_node;
    int edge_count = 0;
    int i;

    auto& device_ctx = g_vpr_ctx.device();
    auto const& segment_inf = device_ctx.arch->Segments;

    VTR_ASSERT(chan_type == CHANX || chan_type == CHANY);
    auto const& chan_detail = chan_details[x_coord][y_coord];
    for (auto wire_type : wire_type_sizes) {
        if (wire_type.second.length == 1)
            continue;

        for (i = 0; i < segment_inf.size(); i++){
            if(wire_type.first == segment_inf[i].name)
                break;
        }
        VTR_ASSERT(i < segment_inf.size());
        std::vector<Wire_Info> real_wire_type;

        if(segment_inf[i].isbend){
            count_bend_wire_type_sizes(chan_detail.data(), wire_type.second.start, wire_type.second.start + wire_type.second.num_wires, real_wire_type);
        }
        else{
            real_wire_type.push_back(wire_type.second);
        }

        for(auto const& wire_inf : real_wire_type){
            int len = wire_inf.length;
            int start = wire_inf.start;
            int num_wires = wire_inf.num_wires;
            int group_size = len * 2; //只考虑单向线
            VTR_ASSERT(num_wires % group_size == 0);
            int group_num = num_wires / group_size;

            for (int ig = 0; ig < group_num; ig++) {
                for (int i_track = start; i_track < start + group_size; i_track++) {
                    auto from_wire = chan_detail[i_track];
                    int actual_len = from_wire.seg_end() - from_wire.seg_start() + 1;
                    if (actual_len == len || from_wire.direction() == to_direction)
                        continue;
                    int pair_len = len - actual_len;
                    for (int i_track2 = start; i_track2 < start + group_size; i_track2++) {
                        if (chan_detail[i_track2].seg_end() - chan_detail[i_track2].seg_start() + 1 == pair_len && chan_detail[i_track2].direction() == to_direction) {
                            from_node = get_rr_node_index(L_rr_node_indices, x_coord, y_coord, chan_type, i_track);
                            to_node = get_rr_node_index(L_rr_node_indices, x_coord, y_coord, chan_type, i_track2);
                            if (from_node == OPEN || to_node == OPEN)
                                break;
                            rr_edges_to_create.emplace_back(from_node, to_node, delayless_switch);
                            ++edge_count;
                            break;
                        }
                    }
                }
                start += group_size;
            }
        }
    }
    return edge_count;
}

int connect_bounds_wires(const t_wire_type_sizes& wire_type_sizes, 
                                const DeviceGrid& grid,
                                const t_chan_details& chan_details_x,
                                const t_chan_details& chan_details_y,
                                const t_rr_node_indices& L_rr_node_indices,
                                t_rr_edge_info_set& rr_edges_to_create,
                                const int delayless_switch) {
    int x_coord, y_coord;
    int from_node, to_node;
    int edge_count = 0;

    //CHANX
    //left bound
    x_coord = 1;
    for (y_coord = 0; y_coord < grid.height() - 1; y_coord++){
        edge_count += build_a_bound_wire(x_coord, y_coord, wire_type_sizes, chan_details_x, L_rr_node_indices, rr_edges_to_create, INC_DIRECTION, CHANX, delayless_switch);
    }

    //right bound
    x_coord = grid.width() -2;
    for (y_coord = 0; y_coord < grid.height() - 1; y_coord++) {
        edge_count += build_a_bound_wire(x_coord, y_coord, wire_type_sizes, chan_details_x, L_rr_node_indices, rr_edges_to_create, DEC_DIRECTION, CHANX, delayless_switch);
    }

    //CHANY
    //bottom bound
    y_coord = 1;
    for (x_coord = 0; x_coord < grid.width() - 1; x_coord++) {
        edge_count += build_a_bound_wire(x_coord, y_coord, wire_type_sizes, chan_details_y, L_rr_node_indices, rr_edges_to_create, INC_DIRECTION, CHANY, delayless_switch);
    }

    //top bound
    y_coord = grid.height() -2;
    for (x_coord = 0; x_coord < grid.width() - 1; x_coord++) {
        edge_count += build_a_bound_wire(x_coord, y_coord, wire_type_sizes, chan_details_y, L_rr_node_indices, rr_edges_to_create, DEC_DIRECTION, CHANY, delayless_switch);
    }

    return edge_count;
}

int build_two_stage_mux_for_gsb(int x_coord, int y_coord, t_rr_edge_info_set& rr_edges_to_create, t_two_stage_mux_map* two_stage_mux_map) {
    int edge_count = 0;

    GSB_Lookup gsb_stage_loc(x_coord, y_coord, e_from_type::NUM_FT_TYPES, NUM_SIDES, NUM_SIDES);
    std::vector<t_stage_mux_inf> first_stages;
    std::vector<t_stage_mux_inf> second_stages;
    auto it = (*two_stage_mux_map).find(gsb_stage_loc);
    if (it != (*two_stage_mux_map).end()) {
        first_stages = (*it).second.first_stages;
        second_stages = (*it).second.second_stages;
    }

    int to_node = OPEN;
    int switch_id;
    for(auto const& first_stage : first_stages){
        to_node = first_stage.to_node;
        if(to_node == OPEN)
            continue;
        switch_id = first_stage.to_switch;
        for(int from_node : first_stage.from_nodes){
            if(from_node == OPEN)
                continue;
            //VTR_LOG("from_node =%d, to_node=%d, sw=%d\n", from_node, to_node, switch_id);
            rr_edges_to_create.emplace_back(from_node, to_node, switch_id);
            ++edge_count;
        }
    }

    for (auto const& second_stage : second_stages) {
        to_node = second_stage.to_node;
        if (to_node == OPEN)
            continue;
        switch_id = second_stage.to_switch;
        for (int from_node : second_stage.from_nodes) {
            if (from_node == OPEN)
                continue;
            //VTR_LOG("from_node =%d, to_node=%d, sw=%d\n", from_node, to_node, switch_id);
            rr_edges_to_create.emplace_back(from_node, to_node, switch_id);
            ++edge_count;
        }
    }

    return edge_count;
}

void filter_midpoint_edges(t_rr_edge_info_set& rr_edges_to_create, const std::vector<t_segment_inf>& segment_inf) {
    auto& device_ctx = g_vpr_ctx.device();
    int from_node, to_node;

    for(auto it = rr_edges_to_create.begin(); it != rr_edges_to_create.end(); ){
        /*auto& from_node = device_ctx.rr_nodes[it->from_node];
        if(from_node.type() == CHANX || from_node.type() == CHANY){
            if (from_node.length() != segment_inf[from_node.seg_id()].length - 1) {
                rr_edges_to_create.erase(it);
                continue;
            }
        }*/

        auto& to_node = device_ctx.rr_nodes[it->to_node];
        if (to_node.type() == CHANX || to_node.type() == CHANY) {
            if (to_node.length() != segment_inf[to_node.seg_id()].length - 1) {
                rr_edges_to_create.erase(it);
                continue;
            }
        }
        it++;
    }
}