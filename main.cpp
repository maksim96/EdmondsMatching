#include <iostream>
#include <queue>
#include "graph.hpp"

using namespace ED;

bool is_outer(NodeId x,std::vector<NodeId>* mu, std::vector<NodeId>* phi){
    return (((*mu)[x] == x) || ((*phi)[(*mu)[x]] != (*mu)[x]));
}

bool is_out_of_forest(NodeId x, std::vector<NodeId>* mu, std::vector<NodeId>* phi){
    return (((*mu)[x] != x) && ((*phi)[x] == x) && ((*phi)[(*mu)[x]] == (*mu)[x]));
}

// In our implementation rho[x] is not necessary the real rho(x)
// instead it is basically like a parent pointer to the next compressed blossom root
// You have to go completely up if you want to have the actual rho value
// The reason for that "lazy" rho is that we don't have to traverse the whole graph at the end of (6) anymore,
// which heavily increases the performance of our algorithm. It's kind of similar to a Union-Find data structure.
NodeId get_rho(std::vector<NodeId> *rho, NodeId x) {
    NodeId xtemp = x;
    while (x != (*rho)[x]) {
        x = (*rho)[x];
    }
    while (xtemp != (*rho)[xtemp]) {
        (*rho)[xtemp] = (*rho)[x];
        xtemp = (*rho)[xtemp];
    }
    return x;
}


//Returns true <=> V(P(x)) cap V(P(y)) != 0
//Additionally, if the paths are not disjoint:
//it sets on_x_r_y_path to contain the P[x,r] and P[y,r] path and finds r
//It assumes that the caller made sure that the on_x_r_y_path array is completely set to false!
//
//If the paths are disjoint, on_x_r_y_path will mark just both paths, but this information won't be used afterwards
bool are_paths_disjoint(NodeId x, NodeId y, NodeId* r, std::vector<NodeId>* mu, std::vector<NodeId>* phi, std::vector<NodeId>* rho, std::vector<bool>* on_x_r_y_path) {
    bool are_disjoint = true;

    //Visit P(x)
    int i = 0;
    (*on_x_r_y_path)[x] = true;
    while ((i % 2 == 0 && (*mu)[x] != x) || (i % 2 == 1 && (*phi)[x] != x)) {
        if (i % 2 == 0) {
            x = (*mu)[x];
        } else {
            x = (*phi)[x];
        }
        (*on_x_r_y_path)[x] = true;
        i++;
    }
    //visit P(y)
    i = 0;
    bool found_r = false;
    while (true) {
        if ((*on_x_r_y_path)[y] && are_disjoint) {
            //found vertex in both paths
            //rest of the path will be the same
            are_disjoint = false;
        }
        (*on_x_r_y_path)[y] = true;
        if (!are_disjoint) {
            //on this place we can actually use the lazy rho[y] instead of get_rho
            //since we just have to now whether rho(y) != y and not compute the real rho(y).
            if (!found_r && (*rho)[y] == y) {
                *r = y;
                found_r = true;
            } else if (found_r) {
                (*on_x_r_y_path)[y] = false;
            }
        }

        //traverse
        if (i % 2 == 0 && y != (*mu)[y]) {
            y = (*mu)[y];
        } else if (i % 2 == 1 && y != (*phi)[y]) {
            y = (*phi)[y];
        } else {
            break;
        }
        i++;

    }

    return are_disjoint;
}

// n = |V(G)|, m =|E(G)|
//prints the graph in dmx format to std out
void print_matching(std::vector<NodeId> *mu, NodeId n) {
    int matching_edges = 0;
    for (NodeId v = 0; v < n; v++) {
        if (v < (*mu)[v]) {
            matching_edges++;
        }
    }
    std::cout << "p edge " << n << " " << matching_edges << std::endl;
    for (NodeId v = 0; v < n; v++) {
        if (v < (*mu)[v]) {
            std::cout << "e " << v + 1 << " " << (*mu)[v] + 1 << std::endl;
        }
    }
}

//Computes a maximum matching of g with the Edmonds' Blossom Algorithm.
//Returns the matching given by the mu vector
std::vector<NodeId> edmonds (Graph g) {
    NodeId n = g.num_nodes();
    //those arrays model exactly mu, phi, scanned(x) and whether a vertex is on P[x,r] union P[y,r]
    std::vector<NodeId> mu(n);
    std::vector<NodeId> phi(n);
    std::vector<bool> scanned(n);
    std::vector<bool> on_x_r_y_path(n);
    //rho array stores the real rho values implicitly, by pointing to the "next" blossom root
    //use get_rho to access the real value
    std::vector<NodeId> rho(n);
    //We use a queue to collect all odd_vertices during the path traversal, so we don't interfere
    //and can update the values afterwards.
    std::queue<NodeId> odd_vertices;
    NodeId r;
    //vertex v and int i for traversal
    NodeId v;
    int i;

    //(1) init
    for (NodeId i = 0; i < n; i++) {
        mu[i] = i;
        phi[i] = i;
        rho[i] = i;
    }

    int matchings_found = 0;

    //Greedy start
    //Just find valid matching edges as long as you can
    for (NodeId i = 0; i < n; i++) {
        if (mu[i] != i) {
            continue;
        }
        for (auto &j: g.node(i).neighbors()) {
            if (mu[j] != j) {
                continue;
            }
            mu[i] = j;
            mu[j] = i;
            matchings_found++;
            break;
        }
    }

    bool found_x = true;
    bool found_y = false;

    while(found_x) {
        //(2)
        found_x = false;
        for (NodeId x = 0; x < n; x++) {
            if (!(is_outer(x, &mu, &phi) && !scanned[x])) {
                continue;
            }
            found_x = true;
            //(3)
            for (auto &y: g.node(x).neighbors()) {
                found_y = false;
                //basically (4)
                if (is_out_of_forest(y, &mu, &phi)){
                    phi[y] = x;
                    found_y = true;
                }
                if (is_outer(y, &mu, &phi) && get_rho(&rho, x) != get_rho(&rho, y)) {
                    found_y = true;
                    if (are_paths_disjoint(x, y, &r, &mu, &phi, &rho, &on_x_r_y_path)) {
                        // (5)
                        //for all v in V(P(x)) union V(P(y)) collect the odd vertices in a queue and clean up the on_x_r_y_path markers directly
                        //traverse P(x)
                        i = 0;
                        v = x;
                        on_x_r_y_path[v] = false;
                        while ((i % 2 == 0 && mu[v] != v) || (i % 2 == 1 && phi[v] != v)) {

                            if (i % 2 == 0) {
                                //traverse
                                v = mu[v];
                                odd_vertices.push(v);
                            } else {
                                //traverse
                                v = phi[v];
                            }
                            i++;
                            on_x_r_y_path[v] = false;
                        }
                        //traverse P(y)
                        v = y;
                        i = 0;
                        on_x_r_y_path[v] = false;
                        while ((i % 2 == 0 && mu[v] != v) || (i % 2 == 1 && phi[v] != v)) {

                            if (i % 2 == 0) {
                                //traverse
                                v = mu[v];
                                odd_vertices.push(v);
                            } else {
                                //traverse
                                v = phi[v];
                            }
                            i++;
                            on_x_r_y_path[v] = false;
                        }

                        //update the odd vertices
                        while (!odd_vertices.empty()) {
                            v = odd_vertices.front();
                            //set new values on odd distance
                            mu[phi[v]] = v;
                            mu[v] = phi[v];
                            odd_vertices.pop();
                        }

                        //do the rest of (5)
                        mu[x] = y;
                        mu[y] = x;

                        //reset scanned, phi, rho
                        for (NodeId v = 0; v < n; v++) {
                            phi[v] = v;
                            rho[v] = v;
                            scanned[v] = false;
                        }
                        break; //aka goto 2
                    } else {
                        // (6)
                        //Overall we travers the paths three times here.
                        //1. To collect the odd vertices
                        //2. Update the rhos only at the paths
                        //3. Clean up the on_x_r_y_path markers
                        //In between 1 and 2 we update the odd vertices
                        //We need this separation as otherwise the updates of mu,phi or rho may interfere the path traversing

                        //for all v in V(P(x)[x,r]) union V(P(y)[y,r]) collect odd vertices
                        //traverse P(x)[x,r]
                        i = 0;
                        v = x;
                        while ((i % 2 == 0 && mu[v] != v && on_x_r_y_path[mu[v]]) || (i % 2 == 1 && phi[v] != v && on_x_r_y_path[phi[v]])) {

                            if (i % 2 == 0) {
                                //traverse
                                v = mu[v];
                                //set new values on odd distance
                                if (get_rho(&rho, phi[v]) != r) {
                                    odd_vertices.push(v);
                                }
                            } else {
                                //traverse
                                v = phi[v];
                            }
                            i++;
                        }

                        //traverse P(y)[y,r]
                        v = y;
                        i = 0;
                        while ((i % 2 == 0 && mu[v] != v && on_x_r_y_path[mu[v]]) || (i % 2 == 1 && phi[v] != v && on_x_r_y_path[phi[v]])) {

                            if (i % 2 == 0) {
                                //traverse
                                v = mu[v];
                                if (get_rho(&rho, phi[v]) != r) {
                                    odd_vertices.push(v);
                                }
                                //set new values on odd distance

                            } else {
                                //traverse
                                v = phi[v];
                            }
                            i++;
                        }

                        //update the odd vertices
                        while (!odd_vertices.empty()) {
                            v = odd_vertices.front();
                            phi[phi[v]] = v;
                            odd_vertices.pop();
                        }

                        //update phi[x], phi[y]
                        if (get_rho(&rho, x) != r) {
                            phi[x] = y;
                        }
                        if (get_rho(&rho, y)  != r) {
                            phi[y] = x;
                        }

                        //Update rho but "lazily"
                        //Update only the rhos on the path and don't traverse the whole graph as in the pseudocode
                        i = 0;
                        v = x;
                        if (on_x_r_y_path[get_rho(&rho, v)]) {
                            rho[v] = r;
                        }
                        while ((i % 2 == 0 && mu[v] != v && on_x_r_y_path[mu[v]]) || (i % 2 == 1 && phi[v] != v && on_x_r_y_path[phi[v]])) {
                            if (i % 2 == 0) {
                                //traverse
                                v = mu[v];
                            } else {
                                //traverse
                                v = phi[v];
                            }
                            i++;
                            if (on_x_r_y_path[get_rho(&rho, v)]) {
                                rho[v] = r;
                            }
                        }

                        v = y;
                        i = 0;
                        if (on_x_r_y_path[get_rho(&rho, v)]) {
                            rho[v] = r;
                        }
                        while ((i % 2 == 0 && mu[v] != v && on_x_r_y_path[mu[v]]) || (i % 2 == 1 && phi[v] != v && on_x_r_y_path[phi[v]])) {

                            if (i % 2 == 0) {
                                //traverse
                                v = mu[v];
                            } else {
                                //traverse
                                v = phi[v];
                            }
                            i++;
                            if (on_x_r_y_path[get_rho(&rho, v)]) {
                                rho[v] = r;
                            }

                        }

                        //Clean up the on_x_r_y_path markers
                        i = 0;
                        v = x;
                        on_x_r_y_path[v] = false;
                        on_x_r_y_path[r] = false;
                        while ((i % 2 == 0 && mu[v] != v && mu[v] != r) || (i % 2 == 1 && phi[v] != v && phi[v] != r)) {

                            if (i % 2 == 0) {
                                //traverse
                                v = mu[v];
                            } else {
                                //traverse
                                v = phi[v];
                            }
                            i++;
                            on_x_r_y_path[v] = false;

                        }



                        v = y;
                        i = 0;
                        on_x_r_y_path[v] = false;
                        while ((i % 2 == 0 && mu[v] != v) || (i % 2 == 1 && phi[v] != v)) {

                            if (i % 2 == 0) {
                                //traverse
                                v = mu[v];
                            } else {
                                //traverse
                                v = phi[v];
                            }
                            i++;
                            on_x_r_y_path[v] = false;

                        }




                    }
                }


            }
            // last part of (3)
            if (!found_y) {
                scanned[x] = true;
                //go to (2) by while loop
            }


        }

    }
    return mu;
}

int main(int argc, char** argv)
{
    if (argc != 2)
    {
        std::cerr << "Wrong number of arguments. Program call: <program_name> <input_graph>" << std::endl;
        return EXIT_FAILURE;
    }

    ED::Graph graph = ED::Graph::build_graph(argv[1]);

    std::vector<NodeId> mu = edmonds(graph);

    print_matching(&mu, graph.num_nodes());
    return EXIT_SUCCESS;







}

