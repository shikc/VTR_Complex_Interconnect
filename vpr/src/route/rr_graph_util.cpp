#include "vtr_memory.h"

#include "vpr_types.h"
#include "vpr_error.h"

#include "globals.h"
#include "rr_graph_util.h"

int seg_index_of_cblock(t_rr_type from_rr_type, int to_node) {
    /* Returns the segment number (distance along the channel) of the connection *
     * box from from_rr_type (CHANX or CHANY) to to_node (IPIN).                 */

    auto& device_ctx = g_vpr_ctx.device();

    if (from_rr_type == CHANX)
        return (device_ctx.rr_nodes[to_node].xlow());
    else
        /* CHANY */
        return (device_ctx.rr_nodes[to_node].ylow());
}

int seg_index_of_sblock(int from_node, int to_node) {
    /* Returns the segment number (distance along the channel) of the switch box *
     * box from from_node (CHANX or CHANY) to to_node (CHANX or CHANY).  The     *
     * switch box on the left side of a CHANX segment at (i,j) has seg_index =   *
     * i-1, while the switch box on the right side of that segment has seg_index *
     * = i.  CHANY stuff works similarly.  Hence the range of values returned is *
     * 0 to device_ctx.grid.width()-1 (if from_node is a CHANX) or 0 to device_ctx.grid.height()-1 (if from_node is a CHANY).   */

    t_rr_type from_rr_type, to_rr_type;

    auto& device_ctx = g_vpr_ctx.device();

    from_rr_type = device_ctx.rr_nodes[from_node].type();
    to_rr_type = device_ctx.rr_nodes[to_node].type();

    if (from_rr_type == CHANX) {
        if (to_rr_type == CHANY) {
            return (device_ctx.rr_nodes[to_node].xlow());
        } else if (to_rr_type == CHANX) {
            if (device_ctx.rr_nodes[to_node].xlow() > device_ctx.rr_nodes[from_node].xlow()) { /* Going right */
                return (device_ctx.rr_nodes[from_node].xhigh());
            } else { /* Going left */
                return (device_ctx.rr_nodes[to_node].xhigh());
            }
        } else {
            vpr_throw(VPR_ERROR_ROUTE, __FILE__, __LINE__,
                      "in seg_index_of_sblock: to_node %d is of type %d.\n",
                      to_node, to_rr_type);
            return OPEN; //Should not reach here once thrown
        }
    }
    /* End from_rr_type is CHANX */
    else if (from_rr_type == CHANY) {
        if (to_rr_type == CHANX) {
            return (device_ctx.rr_nodes[to_node].ylow());
        } else if (to_rr_type == CHANY) {
            if (device_ctx.rr_nodes[to_node].ylow() > device_ctx.rr_nodes[from_node].ylow()) { /* Going up */
                return (device_ctx.rr_nodes[from_node].yhigh());
            } else { /* Going down */
                return (device_ctx.rr_nodes[to_node].yhigh());
            }
        } else {
            vpr_throw(VPR_ERROR_ROUTE, __FILE__, __LINE__,
                      "in seg_index_of_sblock: to_node %d is of type %d.\n",
                      to_node, to_rr_type);
            return OPEN; //Should not reach here once thrown
        }
    }
    /* End from_rr_type is CHANY */
    else {
        vpr_throw(VPR_ERROR_ROUTE, __FILE__, __LINE__,
                  "in seg_index_of_sblock: from_node %d is of type %d.\n",
                  from_node, from_rr_type);
        return OPEN; //Should not reach here once thrown
    }
}

int get_rr_node_index(const t_rr_node_indices& L_rr_node_indices,
                      int x,
                      int y,
                      t_rr_type rr_type,
                      int ptc,
                      e_side side) {
    /*
     * Returns the index of the specified routing resource node.  (x,y) are
     * the location within the FPGA, rr_type specifies the type of resource,
     * and ptc gives the number of this resource.  ptc is the class number,
     * pin number or track number, depending on what type of resource this
     * is.  All ptcs start at 0 and go up to pins_per_clb-1 or the equivalent.
     * There are type->num_class SOURCEs + SINKs, type->num_pins IPINs + OPINs,
     * and max_chan_width CHANX and CHANY (each).
     *
     * Note that for segments (CHANX and CHANY) of length > 1, the segment is
     * given an rr_index based on the (x,y) location at which it starts (i.e.
     * lowest (x,y) location at which this segment exists).
     * This routine also performs error checking to make sure the node in
     * question exists.
     *
     * The 'side' argument only applies to IPIN/OPIN types, and specifies which
     * side of the grid tile the node should be located on. The value is ignored
     * for non-IPIN/OPIN types
     */
    if (rr_type == IPIN || rr_type == OPIN) {
        VTR_ASSERT_MSG(side != NUM_SIDES, "IPIN/OPIN must specify desired side (can not be default NUM_SIDES)");
    } else {
        VTR_ASSERT(rr_type != IPIN && rr_type != OPIN);
        side = SIDES[0];
    }

    int iclass;

    auto& device_ctx = g_vpr_ctx.device();

    VTR_ASSERT(ptc >= 0);
    VTR_ASSERT(x >= 0 && x < int(device_ctx.grid.width()));
    VTR_ASSERT(y >= 0 && y < int(device_ctx.grid.height()));

    t_type_ptr type = device_ctx.grid[x][y].type;

    /* Currently need to swap x and y for CHANX because of chan, seg convention */
    if (CHANX == rr_type) {
        std::swap(x, y);
    }

    /* Start of that block.  */
    const std::vector<int>& lookup = L_rr_node_indices[rr_type][x][y][side];

    /* Check valid ptc num */
    VTR_ASSERT(ptc >= 0);

    switch (rr_type) {
        case SOURCE:
            VTR_ASSERT(ptc < type->num_class);
            VTR_ASSERT(type->class_inf[ptc].type == DRIVER);
            break;

        case SINK:
            VTR_ASSERT(ptc < type->num_class);
            VTR_ASSERT(type->class_inf[ptc].type == RECEIVER);
            break;

        case OPIN:
            /*if(ptc>type->num_pins){
                int a;
                a = ptc;
                int b =type->num_pins;
            }*/
            VTR_ASSERT(ptc < type->num_pins);
            iclass = type->pin_class[ptc];
            /*if(type->class_inf[iclass].type != DRIVER){
                int a;
                a = ptc;
                int b =type->num_pins;
            }*/
            VTR_ASSERT(type->class_inf[iclass].type == DRIVER);
            break;

        case IPIN:
            VTR_ASSERT(ptc < type->num_pins);
            iclass = type->pin_class[ptc];
            VTR_ASSERT(type->class_inf[iclass].type == RECEIVER);
            break;

        case CHANX:
        case CHANY:
            break;

        default:
            vpr_throw(VPR_ERROR_ROUTE, __FILE__, __LINE__,
                      "Bad rr_node passed to get_rr_node_index.\n"
                      "Request for type=%d ptc=%d at (%d, %d).\n",
                      rr_type, ptc, x, y);
    }

    return ((unsigned)ptc < lookup.size() ? lookup[ptc] : -1);
}

int get_medium_rr_node_index(const t_medium_rr_node_indices& L_medium_rr_node_indices,
                             int x,
                             int y,
                             int medium_type,//0-IMUX, 1-OMUX, 2-GSB
                             int num_type_ind,//类型索引，如今默认都只有一种
                             int ptc) {

    VTR_ASSERT( medium_type < 3);
    auto& device_ctx = g_vpr_ctx.device();
    VTR_ASSERT(x >= 0 && x < int(device_ctx.grid.width()));
    VTR_ASSERT(y >= 0 && y < int(device_ctx.grid.height()));
    //t_type_ptr type = device_ctx.grid[x][y].type;
    const std::vector<int>& lookup = L_medium_rr_node_indices[x][y][medium_type][num_type_ind];
    VTR_ASSERT(ptc >= 0);
    return ((unsigned)ptc < lookup.size() ? lookup[ptc] : -1);
}
