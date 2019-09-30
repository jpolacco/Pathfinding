package roadgraph;

import geography.GeographicPoint;
import geography.RoadSegment;

import java.util.ArrayList;
import java.util.List;

/**
 *  The MapNode class wraps both the GeographicPoint of the node and the all the unidirectional out-edges from this
 *  GeographicPoint to other intersections.
 */
public class MapNode {
    private GeographicPoint point;
    private List<MapEdge> edges;

    /**
     * Creates a MapNode object given the intersection and the list of edges.
     * @param gp the intersection of the MapNode
     * @param me the list of unidirectional out-edges for this intersection
     */
    public MapNode(GeographicPoint gp, List<MapEdge> me){
        /*
        Note to myself: Should I make deep copies of these objects here eventually?
        */
        point = gp;
        edges = me;
    }

    /**
     * Adds an edge to this MapNode
     * @param edge  The MapEdge that is being added to this MapNode
     */
    public void addEdge(MapEdge edge){
        edges.add(edge);
    }

    /**
     * Gets the intersection for this MapNode.
     * @return the intersection for this MapNode
     */
    public GeographicPoint getGeographicPoint(){
        return point;
    }

    /**
     * Gets all the out-neighbors for this MapNode
     * @return all the out-neighbors for this MapNode
     */
    public List<GeographicPoint> getNeighbors(){
        List<GeographicPoint> neighbors = new ArrayList<GeographicPoint>();
        for (MapEdge edge : edges){
            neighbors.add(edge.getEndPoint());
        }
        return neighbors;
    }


}
