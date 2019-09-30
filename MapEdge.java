package roadgraph;

import geography.GeographicPoint;


/**
 * The MapEdge class encapsulates both the starting and ending GeographicPoint, as well
 * as the road name and distance.
 */
public class MapEdge {
    private GeographicPoint start;
    private GeographicPoint end;
    private String streetName;
    private double distance;

    /**
     * Creates a MapEdge object given the starting and ending intersections of the edge, street name, and distance.
     * @param start  The starting intersection of the edge
     * @param end  The ending intersection of the edge
     * @param streetName  The name of the street on this edge
     * @param distance  The distance from the start to the end
     */
    public MapEdge(GeographicPoint start, GeographicPoint end, String streetName, double distance ){
        this.start = start;
        this.end = end;
        this.streetName = streetName;
        this.distance = distance;
    }

    /**
     * Get the intersection at the end of this road segment(edge).
     * @return The intersection at the end of this road segment(edge).
     */
    public GeographicPoint getEndPoint(){
        return end;
    }

    /**
     * Get the intersection at the start of this road segment(edge).
     * @return The intersection at the start of this road segment(edge).
     */
    public GeographicPoint getStartPoint(){
        return start;
    }
}
