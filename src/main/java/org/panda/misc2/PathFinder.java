package org.panda.misc2;

import org.pathwaycommons.sif.io.Loader;
import org.pathwaycommons.sif.io.Writer;
import org.pathwaycommons.sif.model.SIFGraph;
import org.pathwaycommons.sif.query.QueryExecutor;
import org.pathwaycommons.sif.util.EdgeAnnotationType;
import org.pathwaycommons.sif.util.EdgeSelector;
import org.pathwaycommons.sif.util.RelationTypeSelector;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;
import java.util.zip.GZIPInputStream;

import static org.pathwaycommons.sif.model.RelationTypeEnum.CONTROLS_EXPRESSION_OF;
import static org.pathwaycommons.sif.model.RelationTypeEnum.CONTROLS_STATE_CHANGE_OF;
import static org.pathwaycommons.sif.util.EdgeAnnotationType.*;

public class PathFinder
{
    public static void main(String[] args) throws IOException {
        // Tell the type and order of edge annotations in the SIF resource
        EdgeAnnotationType[] edgeAnnotTypes = new EdgeAnnotationType[]{
                DATA_SOURCE, PUBMED_IDS, PATHWAY_NAMES, MEDIATORS};

        EdgeAnnotationType[] outTypes = new EdgeAnnotationType[]{MEDIATORS};



        // Initialize loader
        Loader loader = new Loader(edgeAnnotTypes);

        // Load a graph
        SIFGraph graph = loader.load(new GZIPInputStream( new FileInputStream("/Users/ozgun/Downloads/PathwayCommons12.Detailed.hgnc.txt.gz")));

        //--- Perform a query on the graph

        // The query will traverse only two type of relations in this example
        EdgeSelector edgeSelector = new RelationTypeSelector(CONTROLS_STATE_CHANGE_OF, CONTROLS_EXPRESSION_OF);

        // Select a seed
        Set<String> source = new HashSet<>(Arrays.asList("EGF"));
        Set<String> target = new HashSet<>(Arrays.asList("JUN"));

        Set<Object> result = Collections.emptySet();

        for (int limit = 1; result.isEmpty() && limit < 5; limit++)
        {
            // Run the query
            result = QueryExecutor.searchPathsFromTo(graph, edgeSelector, source, target, limit);
        }

        //--- Report results

        // Initialize the writer with the same edge annotation style
        Writer writer = new Writer(false, outTypes);

        // Write results
        writer.write(result, System.out);

    }
}
