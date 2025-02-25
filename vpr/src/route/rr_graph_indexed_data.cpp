#include <cmath> /* Needed only for sqrt call (remove if sqrt removed) */
using namespace std;

#include "vtr_assert.h"
#include "vtr_log.h"
#include "vtr_memory.h"

#include "vpr_types.h"
#include "vpr_error.h"

#include "globals.h"
#include "rr_graph_util.h"
#include "rr_graph2.h"
#include "rr_graph_indexed_data.h"
#include "read_xml_arch_file.h"

/******************* Subroutines local to this module ************************/

static void load_rr_indexed_data_base_costs(int nodes_per_chan,
                                            const t_rr_node_indices& L_rr_node_indices,
                                            enum e_base_cost_type base_cost_type);

static float get_delay_normalization_fac(int nodes_per_chan,
                                         const t_rr_node_indices& L_rr_node_indices);

static void load_rr_indexed_data_T_values(int index_start,
                                          int num_indices_to_load,
                                          t_rr_type rr_type,
                                          int nodes_per_chan,
                                          const t_medium_seg_inf& medium_seg_inf,
                                          bool is_y_medium,
                                          const t_rr_node_indices& L_rr_node_indices,
                                          const t_medium_rr_node_indices& L_medium_rr_node_indices);

static std::vector<size_t> count_rr_segment_types();

/******************** Subroutine definitions *********************************/

/* Allocates the device_ctx.rr_indexed_data array and loads it with appropriate values. *
 * It currently stores the segment type (or OPEN if the index doesn't        *
 * correspond to an CHANX or CHANY type), the base cost of nodes of that     *
 * type, and some info to allow rapid estimates of time to get to a target   *
 * to be computed by the router.                                             *
 *
 * Right now all SOURCES have the same base cost; and similarly there's only *
 * one base cost for each of SINKs, OPINs, and IPINs (four total).  This can *
 * be changed just by allocating more space in the array below and changing  *
 * the cost_index values for these rr_nodes, if you want to make some pins   *
 * etc. more expensive than others.  I give each segment type in an          *
 * x-channel its own cost_index, and each segment type in a y-channel its    *
 * own cost_index.                                                           */
void alloc_and_load_rr_indexed_data(const std::vector<t_segment_inf>& segment_inf,
                                    const t_medium_seg_inf& medium_seg_inf,
                                    const t_rr_node_indices& L_rr_node_indices,
                                    const t_medium_rr_node_indices& L_medium_rr_node_indices,
                                    const int nodes_per_chan,
                                    int wire_to_ipin_switch,
                                    enum e_base_cost_type base_cost_type) {
    int iseg, length, i, index;

    /* reload the num_segment, cut off the bend segment */
    int real_num_segment = 0;
    vector<int> real_length;
    vector<bool> segment_longline;
    vector<int> real_seg_index;
    for (iseg = 0; iseg < segment_inf.size(); iseg++) {
        if (segment_inf[iseg].isbend) {
            int tmp_length = 1;
            for (auto bend_seg = segment_inf[iseg].bend.begin(); bend_seg != segment_inf[iseg].bend.end(); bend_seg++) {
                if (*bend_seg != 0) {
                    real_length.push_back(tmp_length);
                    segment_longline.push_back(segment_inf[iseg].longline);
                    real_seg_index.push_back(iseg);
                    real_num_segment++;
                    tmp_length = 1;
                } else {
                    tmp_length++;
                }
            }
            // the last part
            real_num_segment++;
            real_length.push_back(tmp_length);
            segment_longline.push_back(segment_inf[iseg].longline);
            real_seg_index.push_back(iseg);
        } else {
            real_length.push_back(segment_inf[iseg].length);
            segment_longline.push_back(segment_inf[iseg].longline);
            real_seg_index.push_back(iseg);
            real_num_segment++;
        }
    }

    VTR_ASSERT(real_length.size() == segment_longline.size());
    VTR_ASSERT(real_length.size() == segment_longline.size());

    auto& device_ctx = g_vpr_ctx.mutable_device();
    int num_segment = segment_inf.size();
    int num_rr_indexed_data = CHANX_COST_INDEX_START + (2 * real_num_segment);
    device_ctx.rr_indexed_data.resize(num_rr_indexed_data);

    /* For rr_types that aren't CHANX or CHANY, base_cost is valid, but most     *
     * * other fields are invalid.  For IPINs, the T_linear field is also valid;   *
     * * all other fields are invalid.  For SOURCES, SINKs and OPINs, all fields   *
     * * other than base_cost are invalid. Mark invalid fields as OPEN for safety. */

    for (i = SOURCE_COST_INDEX; i <= IPIN_COST_INDEX; i++) {
        device_ctx.rr_indexed_data[i].ortho_cost_index = OPEN;
        device_ctx.rr_indexed_data[i].seg_index = OPEN;
        device_ctx.rr_indexed_data[i].inv_length = OPEN;
        device_ctx.rr_indexed_data[i].T_linear = OPEN;
        device_ctx.rr_indexed_data[i].T_quadratic = OPEN;
        device_ctx.rr_indexed_data[i].C_load = OPEN;
    }
    device_ctx.rr_indexed_data[IPIN_COST_INDEX].T_linear = device_ctx.rr_switch_inf[wire_to_ipin_switch].Tdel;

    /* X-directed segments. */
    for (iseg = 0; iseg < real_num_segment; iseg++) {
        index = CHANX_COST_INDEX_START + iseg;

        if ((index + real_num_segment) >= (int)device_ctx.rr_indexed_data.size()) {
            device_ctx.rr_indexed_data[index].ortho_cost_index = index;
        } else {
            device_ctx.rr_indexed_data[index].ortho_cost_index = index + real_num_segment;
        }

        if (segment_inf[iseg].longline)
            length = device_ctx.grid.width();
        else
            length = std::min<int>(real_length[iseg], device_ctx.grid.width());

        device_ctx.rr_indexed_data[index].inv_length = 1. / length;
        device_ctx.rr_indexed_data[index].seg_index = real_seg_index[iseg];
    }
    load_rr_indexed_data_T_values(CHANX_COST_INDEX_START, real_num_segment, CHANX,
                                  nodes_per_chan, medium_seg_inf, false, L_rr_node_indices, L_medium_rr_node_indices);

    /* Y-directed segments. */
    for (iseg = 0; iseg < real_num_segment; iseg++) {
        index = CHANX_COST_INDEX_START + real_num_segment + iseg;

        if ((index - real_num_segment) < CHANX_COST_INDEX_START) {
            device_ctx.rr_indexed_data[index].ortho_cost_index = index;
        } else {
            device_ctx.rr_indexed_data[index].ortho_cost_index = index - real_num_segment;
        }

        if (segment_inf[iseg].longline)
            length = device_ctx.grid.height();
        else
            length = std::min<int>(real_length[iseg], device_ctx.grid.height());

        device_ctx.rr_indexed_data[index].inv_length = 1. / length;
        device_ctx.rr_indexed_data[index].seg_index = real_seg_index[iseg];
    }
    load_rr_indexed_data_T_values((CHANX_COST_INDEX_START + real_num_segment), real_num_segment, CHANY,
                                    nodes_per_chan, medium_seg_inf, true, L_rr_node_indices, L_medium_rr_node_indices);

    load_rr_indexed_data_base_costs(nodes_per_chan, L_rr_node_indices,
                                    base_cost_type);
}

void load_rr_index_segments(const int num_segment, const std::vector<t_segment_inf>& segment_inf) {
    auto& device_ctx = g_vpr_ctx.mutable_device();
    int iseg, i, index;

    for (i = SOURCE_COST_INDEX; i <= IPIN_COST_INDEX; i++) {
        device_ctx.rr_indexed_data[i].seg_index = OPEN;
    }

    /* X-directed segments. */
    int index_offset = 0;
    for (iseg = 0; iseg < num_segment; iseg++) {
        index = CHANX_COST_INDEX_START + iseg + index_offset;
        /* get the real index */
        if (segment_inf[iseg].isbend) {
            device_ctx.rr_indexed_data[index].seg_index = iseg;
            for (auto bend_seg = segment_inf[iseg].bend.begin(); bend_seg != segment_inf[iseg].bend.end(); bend_seg++) {
                if (*bend_seg != 0) {
                    index++;
                    index_offset++;
                    device_ctx.rr_indexed_data[index].seg_index = iseg;
                }
            }
        } else {
            device_ctx.rr_indexed_data[index].seg_index = iseg;
        }
    }

    /* Y-directed segments. */
    for (iseg = 0; iseg < num_segment; iseg++) {
        index = CHANX_COST_INDEX_START + num_segment + index_offset + iseg;
        /* get the real index */
        if (segment_inf[iseg].isbend) {
            device_ctx.rr_indexed_data[index].seg_index = iseg;
            for (auto bend_seg = segment_inf[iseg].bend.begin(); bend_seg != segment_inf[iseg].bend.end(); bend_seg++) {
                if (*bend_seg != 0) {
                    index++;
                    index_offset++;
                    device_ctx.rr_indexed_data[index].seg_index = iseg;
                }
            }
        } else {
            device_ctx.rr_indexed_data[index].seg_index = iseg;
        }
    }
}

static void load_rr_indexed_data_base_costs(int nodes_per_chan,
                                            const t_rr_node_indices& L_rr_node_indices,
                                            enum e_base_cost_type base_cost_type) {
    /* Loads the base_cost member of device_ctx.rr_indexed_data according to the specified *
     * base_cost_type.                                                          */

    float delay_normalization_fac;
    size_t index;

    auto& device_ctx = g_vpr_ctx.mutable_device();

    if (base_cost_type == DEMAND_ONLY || base_cost_type == DEMAND_ONLY_NORMALIZED_LENGTH) {
        delay_normalization_fac = 1.;
    } else {
        delay_normalization_fac = get_delay_normalization_fac(nodes_per_chan, L_rr_node_indices);
    }

    device_ctx.rr_indexed_data[SOURCE_COST_INDEX].base_cost = delay_normalization_fac;
    device_ctx.rr_indexed_data[SINK_COST_INDEX].base_cost = 0.;
    device_ctx.rr_indexed_data[OPIN_COST_INDEX].base_cost = delay_normalization_fac;
    device_ctx.rr_indexed_data[IPIN_COST_INDEX].base_cost = 0.95 * delay_normalization_fac;

    auto rr_segment_counts = count_rr_segment_types();
    size_t total_segments = std::accumulate(rr_segment_counts.begin(), rr_segment_counts.end(), 0u);

    /* Load base costs for CHANX and CHANY segments */

    //Future Work: Since we can now have wire types which don't connect to IPINs,
    //             perhaps consider lowering cost of wires which connect to IPINs
    //             so they get explored earlier (same rational as lowering IPIN costs)

    for (index = CHANX_COST_INDEX_START; index < device_ctx.rr_indexed_data.size(); index++) {
        if (base_cost_type == DELAY_NORMALIZED || base_cost_type == DEMAND_ONLY) {
            device_ctx.rr_indexed_data[index].base_cost = delay_normalization_fac;

        } else if (base_cost_type == DELAY_NORMALIZED_LENGTH || base_cost_type == DEMAND_ONLY_NORMALIZED_LENGTH) {
            device_ctx.rr_indexed_data[index].base_cost = delay_normalization_fac / device_ctx.rr_indexed_data[index].inv_length;

        } else if (base_cost_type == DELAY_NORMALIZED_FREQUENCY) {
            int seg_index = device_ctx.rr_indexed_data[index].seg_index;
            float freq_fac = float(rr_segment_counts[seg_index]) / total_segments;

            device_ctx.rr_indexed_data[index].base_cost = delay_normalization_fac / freq_fac;

        } else if (base_cost_type == DELAY_NORMALIZED_LENGTH_FREQUENCY) {
            int seg_index = device_ctx.rr_indexed_data[index].seg_index;
            float freq_fac = float(rr_segment_counts[seg_index]) / total_segments;

            //Base cost = delay_norm / (len * freq)
            //device_ctx.rr_indexed_data[index].base_cost = delay_normalization_fac / ((1. / device_ctx.rr_indexed_data[index].inv_length) * freq_fac);

            //Base cost = (delay_norm * len) * (1 + (1-freq))
            device_ctx.rr_indexed_data[index].base_cost = (delay_normalization_fac / device_ctx.rr_indexed_data[index].inv_length) * (1 + (1 - freq_fac));

        } else {
            VPR_THROW(VPR_ERROR_ROUTE, "Unrecognized base cost type");
        }
    }

    /* Save a copy of the base costs -- if dynamic costing is used by the     *
     * router, the base_cost values will get changed all the time and being   *
     * able to restore them from a saved version is useful.                   */

    for (index = 0; index < device_ctx.rr_indexed_data.size(); index++) {
        device_ctx.rr_indexed_data[index].saved_base_cost = device_ctx.rr_indexed_data[index].base_cost;
    }
}

static std::vector<size_t> count_rr_segment_types() {
    std::vector<size_t> rr_segment_type_counts;

    auto& device_ctx = g_vpr_ctx.device();

    for (size_t inode = 0; inode < device_ctx.rr_nodes.size(); ++inode) {
        if (device_ctx.rr_nodes[inode].type() != CHANX && device_ctx.rr_nodes[inode].type() != CHANY) continue;

        int cost_index = device_ctx.rr_nodes[inode].cost_index();

        int seg_index = device_ctx.rr_indexed_data[cost_index].seg_index;

        VTR_ASSERT(seg_index != OPEN);

        if (seg_index >= int(rr_segment_type_counts.size())) {
            rr_segment_type_counts.resize(seg_index + 1, 0);
        }
        VTR_ASSERT(seg_index < int(rr_segment_type_counts.size()));

        ++rr_segment_type_counts[seg_index];
    }

    return rr_segment_type_counts;
}

static float get_delay_normalization_fac(int nodes_per_chan,
                                         const t_rr_node_indices& L_rr_node_indices) {
    /* Returns the average delay to go 1 CLB distance along a wire.  */

    const int clb_dist = 3; /* Number of CLBs I think the average conn. goes. */

    int inode, itrack, cost_index;
    float Tdel, Tdel_sum, frac_num_seg;

    auto& device_ctx = g_vpr_ctx.device();

    Tdel_sum = 0.;

    for (itrack = 0; itrack < nodes_per_chan; itrack++) {
        inode = find_average_rr_node_index(device_ctx.grid.width(), device_ctx.grid.height(), CHANX, itrack,
                                           L_rr_node_indices);
        if (inode == -1)
            continue;
        cost_index = device_ctx.rr_nodes[inode].cost_index();
        frac_num_seg = clb_dist * device_ctx.rr_indexed_data[cost_index].inv_length;
        Tdel = frac_num_seg * device_ctx.rr_indexed_data[cost_index].T_linear
               + frac_num_seg * frac_num_seg
                     * device_ctx.rr_indexed_data[cost_index].T_quadratic;
        Tdel_sum += Tdel / (float)clb_dist;
    }

    for (itrack = 0; itrack < nodes_per_chan; itrack++) {
        inode = find_average_rr_node_index(device_ctx.grid.width(), device_ctx.grid.height(), CHANY, itrack,
                                           L_rr_node_indices);
        if (inode == -1)
            continue;
        cost_index = device_ctx.rr_nodes[inode].cost_index();
        frac_num_seg = clb_dist * device_ctx.rr_indexed_data[cost_index].inv_length;
        Tdel = frac_num_seg * device_ctx.rr_indexed_data[cost_index].T_linear
               + frac_num_seg * frac_num_seg
                     * device_ctx.rr_indexed_data[cost_index].T_quadratic;
        Tdel_sum += Tdel / (float)clb_dist;
    }

    return (Tdel_sum / (2. * nodes_per_chan));
}

static void load_rr_indexed_data_T_values(int index_start,
                                          int num_indices_to_load,
                                          t_rr_type rr_type,
                                          int nodes_per_chan,
                                          const t_medium_seg_inf& medium_seg_inf,
                                          bool is_y_medium,
                                          const t_rr_node_indices& L_rr_node_indices,
                                          const t_medium_rr_node_indices& L_medium_rr_node_indices) {
    /* Loads the average propagation times through segments of each index type  *
     * for either all CHANX segment types or all CHANY segment types.  It does  *
     * this by looking at all the segments in one channel in the middle of the  *
     * array and averaging the R and C values of all segments of the same type  *
     * and using them to compute average delay values for this type of segment. */

    int itrack, inode, cost_index;
    float *C_total, *R_total;                /* [0..device_ctx.rr_indexed_data.size() - 1] */
    double *switch_R_total, *switch_T_total; /* [0..device_ctx.rr_indexed_data.size() - 1] */
    short* switches_buffered;
    int* num_nodes_of_index; /* [0..device_ctx.rr_indexed_data.size() - 1] */
    float Rnode, Cnode, Rsw, Tsw;

    auto& device_ctx = g_vpr_ctx.mutable_device();

    num_nodes_of_index = (int*)vtr::calloc(device_ctx.rr_indexed_data.size(), sizeof(int));
    C_total = (float*)vtr::calloc(device_ctx.rr_indexed_data.size(), sizeof(float));
    R_total = (float*)vtr::calloc(device_ctx.rr_indexed_data.size(), sizeof(float));

    /* August 2014: Not all wire-to-wire switches connecting from some wire segment will
     * necessarily have the same delay. i.e. a mux with less inputs will have smaller delay
     * than a mux with a greater number of inputs. So to account for these differences we will
     * get the average R/Tdel values by first averaging them for a single wire segment (first
     * for loop below), and then by averaging this value over all wire segments in the channel
     * (second for loop below) */
    switch_R_total = (double*)vtr::calloc(device_ctx.rr_indexed_data.size(), sizeof(double));
    switch_T_total = (double*)vtr::calloc(device_ctx.rr_indexed_data.size(), sizeof(double));
    switches_buffered = (short*)vtr::calloc(device_ctx.rr_indexed_data.size(), sizeof(short));

    /* initialize switches_buffered array */
    for (int i = index_start; i < index_start + num_indices_to_load; i++) {
        switches_buffered[i] = UNDEFINED;
    }

    /* Get average C and R values for all the segments of this type in one      *
     * channel segment, near the middle of the fpga.                            */
    int medium_num = medium_seg_inf.gsb_mdeium_num_xory;

    for (itrack = 0; itrack < nodes_per_chan + medium_num; itrack++) {
        if(itrack < nodes_per_chan){
            inode = find_average_rr_node_index(device_ctx.grid.width(), device_ctx.grid.height(), rr_type, itrack,
                                           L_rr_node_indices);
        }
        else{
            int medium_index = itrack - nodes_per_chan;
            int medium_type = 0;//0-imux, 1-omux, 2-gsb
            if (medium_index < medium_seg_inf.gsb_mdeium_num_xory){
                medium_index = is_y_medium ? medium_index + medium_seg_inf.gsb_mdeium_num_xory : medium_index;
            }




            inode = find_average_medium_rr_node_index(device_ctx.grid.width(), device_ctx.grid.height(), medium_type, medium_index, L_medium_rr_node_indices);
        }
        if (inode == -1)
            continue;
        cost_index = device_ctx.rr_nodes[inode].cost_index();
        num_nodes_of_index[cost_index]++;
        float cc = device_ctx.rr_nodes[inode].C();
        float ct = C_total[cost_index];
        C_total[cost_index] += device_ctx.rr_nodes[inode].C();
        R_total[cost_index] += device_ctx.rr_nodes[inode].R();

        /* get average switch parameters */
        int num_edges = device_ctx.rr_nodes[inode].num_edges();
        double avg_switch_R = 0;
        double avg_switch_T = 0;
        int num_switches = 0;
        short buffered = UNDEFINED;
        for (int iedge = 0; iedge < num_edges; iedge++) {
            int to_node_index = device_ctx.rr_nodes[inode].edge_sink_node(iedge);
            /* want to get C/R/Tdel of switches that connect this track segment to other track segments */
            //TODO: Test if we should count in the bend segment (xibo)
            if (device_ctx.rr_nodes[to_node_index].type() == CHANX || device_ctx.rr_nodes[to_node_index].type() == CHANY) {
                int switch_index = device_ctx.rr_nodes[inode].edge_switch(iedge);
                avg_switch_R += device_ctx.rr_switch_inf[switch_index].R;
                avg_switch_T += device_ctx.rr_switch_inf[switch_index].Tdel;

                num_switches++;
            }
        }

        if (num_switches == 0) {
            VTR_LOG_WARN("Track %d had no switches\n", itrack);
            continue;
        }
        VTR_ASSERT(num_switches > 0);

        avg_switch_R /= num_switches;
        avg_switch_T /= num_switches;
        switch_R_total[cost_index] += avg_switch_R;
        switch_T_total[cost_index] += avg_switch_T;

        if (buffered == UNDEFINED) {
            /* this segment does not have any outgoing edges to other general routing wires */
            continue;
        }

        /* need to make sure all wire switches of a given wire segment type have the same 'buffered' value */
        if (switches_buffered[cost_index] == UNDEFINED) {
            switches_buffered[cost_index] = buffered;
        } else {
            if (switches_buffered[cost_index] != buffered) {
                vpr_throw(VPR_ERROR_ARCH, __FILE__, __LINE__,
                          "Expecting all wire-to-wire switches of wire segments with cost index (%d) to have same 'buffered' value (%d), but found segment switch with different 'buffered' value (%d)\n", cost_index, switches_buffered[cost_index], buffered);
            }
        }
    }

    for (cost_index = index_start;
         cost_index < index_start + num_indices_to_load; cost_index++) {
        if (num_nodes_of_index[cost_index] == 0) { /* Segments don't exist. */
            device_ctx.rr_indexed_data[cost_index].T_linear = OPEN;
            device_ctx.rr_indexed_data[cost_index].T_quadratic = OPEN;
            device_ctx.rr_indexed_data[cost_index].C_load = OPEN;
        } else {
            Rnode = R_total[cost_index] / num_nodes_of_index[cost_index];
            Cnode = C_total[cost_index] / num_nodes_of_index[cost_index];
            Rsw = (float)switch_R_total[cost_index] / num_nodes_of_index[cost_index];
            Tsw = (float)switch_T_total[cost_index] / num_nodes_of_index[cost_index];

            if (switches_buffered[cost_index]) {
                device_ctx.rr_indexed_data[cost_index].T_linear = Tsw + Rsw * Cnode
                                                                  + 0.5 * Rnode * Cnode;
                device_ctx.rr_indexed_data[cost_index].T_quadratic = 0.;
                device_ctx.rr_indexed_data[cost_index].C_load = 0.;
            } else { /* Pass transistor */
                device_ctx.rr_indexed_data[cost_index].C_load = Cnode;

                /* See Dec. 23, 1997 notes for deriviation of formulae. */

                device_ctx.rr_indexed_data[cost_index].T_linear = Tsw + 0.5 * Rsw * Cnode;
                device_ctx.rr_indexed_data[cost_index].T_quadratic = (Rsw + Rnode) * 0.5
                                                                     * Cnode;
            }
        }
    }

    free(num_nodes_of_index);
    free(C_total);
    free(R_total);
    free(switch_R_total);
    free(switch_T_total);
    free(switches_buffered);
}
