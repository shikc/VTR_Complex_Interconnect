#ifndef BUILD_SWITCHBLOCKS_H
#define BUILD_SWITCHBLOCKS_H

#include <unordered_map>
#include <vector>
#include <random>
#include "physical_types.h"
#include "vpr_types.h"
#include "device_grid.h"
#include "globals.h"
#include "read_xml_arch_file.h"
#include "vtr_random.h"
#include "rr_graph_util.h"


const t_chan_details& index_into_correct_chan(int tile_x, int tile_y, enum e_side side, const t_chan_details& chan_details_x, const t_chan_details& chan_details_y, int* chan_x, int* chan_y, t_rr_type* chan_type);

/* checks whether the specified coordinates are out of bounds */
bool coords_out_of_bounds(const DeviceGrid& grid, int x_coord, int y_coord, e_rr_type chan_type);

t_type_ptr find_pb_type(std::string pb_type_name);

class Wire_Info {
  public:
    int length;    /* the length of this type of wire segment in tiles */
    int num_wires; /* total number of wires in a channel segment (basically W) */
    int start;     /* the wire index at which this type starts in the channel segment (0..W-1) */

    void set(int len, int wires, int st) {
        length = len;
        num_wires = wires;
        start = st;
    }
    Wire_Info() {
        this->set(0, 0, 0);
    }
    Wire_Info(int len, int wires, int st) {
        this->set(len, wires, st);
    }
};

/************ Classes ************/
/* contains info about a wire segment type */

struct t_wire_switchpoint {
    int wire;        //Wire index within the channel
    int switchpoint; //Switchpoint of the wire
};

/************ Typedefs ************/
/* Used to get info about a given wire type based on the name */
typedef std::map<std::string, Wire_Info> t_wire_type_sizes;

/* Counts the number of wires in each wire type in the specified channel */
void count_wire_type_sizes(const t_chan_seg_details* channel, int nodes_per_chan, t_wire_type_sizes* wire_type_sizes);
/************ Classes, structs, typedefs ************/

/* Holds the coordinates of a switch block source connection. Used to index into a
 * map which specifies which destination wire segments this source wire should		//TODO: what data structure does this index to?
 * connect to */
class Switchblock_Lookup {
  public:
    int x_coord; /* x coordinate of switchblock connection */ //TODO: redundant comment?? add range
    int y_coord;                                              /* y coordinate of switchblock connection */
    e_side from_side;                                         /* source side of switchblock connection */
    e_side to_side;                                           /* destination side of switchblock connection */

    /* Empty constructor initializes everything to 0 */
    Switchblock_Lookup() {
        x_coord = y_coord = -1; //TODO: use set function
    }

    /* Constructor for initializing member variables */
    Switchblock_Lookup(int set_x, int set_y, e_side set_from, e_side set_to) {
        this->set_coords(set_x, set_y, set_from, set_to); //TODO: use set function
    }

    /* Function for setting the segment coordinates */
    void set_coords(int set_x, int set_y, e_side set_from, e_side set_to) {
        x_coord = set_x;
        y_coord = set_y;
        from_side = set_from;
        to_side = set_to;
    }

    /* Overload == operator which is used by std::unordered_map */
    bool operator==(const Switchblock_Lookup& obj) const {
        bool result;
        if (x_coord == obj.x_coord && y_coord == obj.y_coord
            && from_side == obj.from_side && to_side == obj.to_side) {
            result = true;
        } else {
            result = false;
        }
        return result;
    }
};

struct t_hash_Switchblock_Lookup {
    size_t operator()(const Switchblock_Lookup& obj) const noexcept {
        //TODO: use vtr::hash_combine
        size_t result;
        result = ((((std::hash<int>()(obj.x_coord)
                     ^ std::hash<int>()(obj.y_coord) << 10)
                    ^ std::hash<int>()((int)obj.from_side) << 20)
                   ^ std::hash<int>()((int)obj.to_side) << 30));
        return result;
    }
};

/* contains the index of the destination wire segment within a channel
 * and the index of the switch used to connect to it */
struct t_switchblock_edge {
    short from_wire;
    short to_wire;
    short switch_ind;
};

/* Switchblock connections are made as [x][y][from_side][to_side][from_wire_ind].
 * The Switchblock_Lookup class specifies these dimensions.
 * Furthermore, a source_wire at a given 5-d coordinate may connect to multiple destination wires so the value
 * of the map is a vector of destination wires.
 * A matrix specifying connections for all switchblocks in an FPGA would be sparse and possibly very large
 * so we use an unordered map to take advantage of the sparsity. */
typedef std::unordered_map<Switchblock_Lookup, std::vector<t_switchblock_edge>, t_hash_Switchblock_Lookup> t_sb_connection_map;

/************ Functions ************/

/* allocate and build switch block permutation map */
t_sb_connection_map* alloc_and_load_switchblock_permutations(const t_chan_details& chan_details_x, const t_chan_details& chan_details_y, const DeviceGrid& grid, std::vector<t_switchblock_inf> switchblocks, t_chan_width* nodes_per_chan, enum e_directionality directionality, vtr::RandState& rand_state);

/* deallocates switch block connections sparse array */
void free_switchblock_permutations(t_sb_connection_map* sb_conns);


//for GSB arch
enum class e_model_type{
    M_IMUX = 0,
    M_OMUX,
    M_GSB,
    NUM_M_TYPES
};

class GSB_Lookup : public Switchblock_Lookup{
  public:
    e_from_type from_type;                                          

    GSB_Lookup() {
        x_coord = y_coord = -1; 
    }

    GSB_Lookup(int set_x, int set_y, e_from_type set_from_type, e_side set_from, e_side set_to) {
        this->set_coords(set_x, set_y, set_from_type, set_from, set_to); //TODO: use set function
    }

    void set_coords(int set_x, int set_y, e_from_type set_from_type, e_side set_from, e_side set_to) {
        x_coord = set_x;
        y_coord = set_y;
        from_type = set_from_type;
        from_side = set_from;
        to_side = set_to;
    }

    /* Overload == operator which is used by std::unordered_map */
    bool operator==(const GSB_Lookup& obj) const {
        bool result;
        if (x_coord == obj.x_coord && y_coord == obj.y_coord && from_type == obj.from_type
            && from_side == obj.from_side && to_side == obj.to_side) {
            result = true;
        } else {
            result = false;
        }
        return result;
    }
};

struct t_hash_GSB_Lookup {
    size_t operator()(const GSB_Lookup& obj) const noexcept {
        size_t result;
        result = (((((std::hash<int>()(obj.x_coord)
                     ^ std::hash<int>()(obj.y_coord) << 10)
                     ^ std::hash<int>()((int)obj.from_side) << 20)
                    ^ std::hash<int>()((int)obj.to_side) << 30)
                   ^ std::hash<int>()((int)obj.from_type) << 40));
        return result;
    }
};

struct t_gsb_edge {
    short from_wireOrpin;
    short to_wireOrpin;
    short switch_ind;
};

typedef std::unordered_map<GSB_Lookup, std::vector<t_gsb_edge>, t_hash_GSB_Lookup> t_gsb_connection_map;

//for gsb multistage mux***************************************
struct t_medium_node_inf{
    int to_index;
    int medium_node;
    std::vector<int> medium_from_nodes;
};

struct t_multistage_mux_map_inf{
    t_medium_node_inf reuse_medium_node;

    std::vector<int> to_nodes;
    std::vector<int> to_switches;
    std::vector<t_medium_node_inf> from_normal_muxes;
};

struct t_multistage_mux_inf_gsbloc {
    std::unordered_map<int, std::pair<int, std::vector<int>>> to_node_loc_map; //to_node, {index of multistage_mux_map_infs, index of from_normal_muxes}
    std::vector<t_multistage_mux_map_inf> multistage_mux_map_infs;
};

typedef std::unordered_map<GSB_Lookup, t_multistage_mux_inf_gsbloc, t_hash_GSB_Lookup> t_gsb_multistage_mux_map;

/************ Functions ************/
t_gsb_connection_map* alloc_and_load_gsb_permutations(const t_chan_details& chan_details_x,
                                                      const t_chan_details& chan_details_y,
                                                      const DeviceGrid& grid,
                                                      std::vector<t_gsb_inf> gsb_fs,
                                                      t_chan_width* nodes_per_chan,
                                                      t_gsb_multistage_mux_map* gsb_multistage_mux_map,
                                                      std::map<std::string, std::vector<std::map<int, int>>> other_pbpin_to_clbpin_map);
void free_gsb_permutations(t_gsb_connection_map * conns);

struct GSB_Side{
    e_side from_side;
    e_side to_side;

    bool operator==(const GSB_Side& obj) const {
        bool result;
        if (from_side == obj.from_side && to_side == obj.to_side) {
            result = true;
        } else {
            result = false;
        }
        return result;
    }
};

struct t_hash_GSB_side {
    size_t operator()(const GSB_Side& obj) const noexcept {
        size_t result;
        result = ((std::hash<int>()(obj.from_side)
                     ^ std::hash<int>()(obj.to_side) << 10));
        return result;
    }
};

struct GSB_detail_conn{
    int from_output;
    int to_input;
};

struct GSB_detail_conn_pb : public GSB_detail_conn{
    std::string port_name;
};

//for gsb twostage mux **************************************************************************
struct t_stage_mux_inf{
    int to_node;
    std::vector<int> from_nodes;
    int to_switch;

    //int gsb_medium_id = 0;
};

struct t_two_stage_mux_inf{
    std::vector<t_stage_mux_inf> first_stages;
    std::vector<t_stage_mux_inf> second_stages;
};

typedef std::unordered_map<GSB_Lookup, t_two_stage_mux_inf, t_hash_GSB_Lookup> t_two_stage_mux_map;

struct t_hash_potential_wire_Lookup {
    size_t operator()(const std::pair<e_side, t_wire_switchpoints>& obj) const noexcept {
        size_t result;
        result = ((((std::hash<int>()(obj.first)
                     ^ std::hash<std::string>()(obj.second.segment_name) << 10)
                    ^ std::hash<int>()((int)*(obj.second.switchpoints.begin() - 1)) << 20)
                   ^ std::hash<int>()((int)*(obj.second.switchpoints.end() - 1)) << 30));
        return result;
    }
};

struct t_potential_wires {
    std::vector<t_wire_switchpoint> potential_wires;
    int x;
    int y;
};
typedef std::unordered_map<std::pair<e_side, t_wire_switchpoints>, std::vector<t_potential_wires>, t_hash_potential_wire_Lookup> t_potential_wires_map;

t_two_stage_mux_map* alloc_and_load_two_stage_map(const t_chan_details& chan_details_x,
                                                                                   const t_chan_details& chan_details_y,
                                                                                   const DeviceGrid& grid,
                                                                                   std::vector<t_gsb_inf> gsb_fs,
                                                                                   const std::vector<t_segment_inf>& segment_inf,
                                                                                   t_chan_width* nodes_per_chan,
                                                                                   const std::map<std::string, std::vector<std::map<int, int>>>& other_pbpin_to_clbpin_map,
                                                                                   const std::map<std::string, int>& gsb_stage_mux_medium_map,
                                                                                   const std::map<int, int>& clbpin_medium_map,
                                                                                   int wire_to_arch_ipin_switch);

//for imux arch ***********************************************
t_gsb_connection_map* alloc_and_load_imux_permutations(const t_chan_details& chan_details_x, const t_chan_details& chan_details_y, const DeviceGrid& grid, std::vector<t_imux_inf> imux_fs, t_chan_width* nodes_per_chan, int wire_to_arch_ipin_switch);

//for omux arch ***********************************************
struct t_omux_edge {
    short from_pb_pin;
    short to_omux_output;
    short switch_ind;
};

typedef std::unordered_map<std::string, std::vector<std::vector<t_omux_edge>>> t_omux_connection_map;
 
t_omux_connection_map* alloc_and_load_omux_permutations(std::vector<t_omux_inf> omux_fs);
void free_omux_permutations(t_omux_connection_map * conns);

template <typename T>
void free_permutations(T* conns) {
    conns->clear();
    delete conns;
    conns = nullptr;
    /* the switch block unordered_map can get quite large and it doesn't seem like the program
     * is interested in releasing the memory back to the OS after the map is cleared.
     * calling malloc_trim forces the program to give unused heap space back to the OS.
     * this significantly reduces memory usage during the routing stage when running multiple
     * large benchmark circuits in parallel. */
    //vtr::malloc_trim(0);
    return;
}

constexpr bool verbose = false;
constexpr bool verbose_imux = false;
#endif
