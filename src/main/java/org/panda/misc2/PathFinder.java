package org.panda.misc2;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.pathwaycommons.sif.io.Loader;
import org.pathwaycommons.sif.model.CustomRelationType;
import org.pathwaycommons.sif.model.SIFEdge;
import org.pathwaycommons.sif.model.SIFGraph;
import org.pathwaycommons.sif.query.Direction;
import org.pathwaycommons.sif.query.QueryExecutor;
import org.pathwaycommons.sif.util.EdgeAnnotationType;
import org.pathwaycommons.sif.util.EdgeSelector;
import org.pathwaycommons.sif.util.RelationTypeSelector;

import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;

import static org.pathwaycommons.sif.model.RelationTypeEnum.INTERACTS_WITH;
import static org.pathwaycommons.sif.model.RelationTypeEnum.IN_COMPLEX_WITH;
import static org.pathwaycommons.sif.model.SignedTypeEnum.*;
import static org.pathwaycommons.sif.util.EdgeAnnotationType.*;

public class PathFinder
{
    public static void main(String[] args) throws IOException {
        // Tell the type and order of edge annotations in the SIF resource
        EdgeAnnotationType[] cpAnnotTypes = new EdgeAnnotationType[]{MEDIATORS, SITES};
        EdgeAnnotationType[] ppiAnnotTypes = new EdgeAnnotationType[]{};

        // Initialize loader
        Loader cpLoader = new Loader(cpAnnotTypes);
        Loader ppiLoader = new Loader(ppiAnnotTypes);

        // Load graphs
        SIFGraph cpGraph = cpLoader.load(new FileInputStream("/home/ozgunbabur/Documents/causal-priors.txt"));
        SIFGraph klGraph = cpLoader.load(new FileInputStream("/home/ozgunbabur/Data/KinaseLibrary/kinase-library.sif"));

        klGraph.getAllEdges().forEach(cpGraph::add);

        SIFGraph ppiGraph = ppiLoader.load(new FileInputStream("/home/ozgunbabur/Documents/PathwayCommons12.All.hgnc.sif"));

        //--- Perform a query on the graph

        // The query will traverse only two type of relations in this example
        EdgeSelector cpEdgeSelector = new RelationTypeSelector(PHOSPHORYLATES, DEPHOSPHORYLATES, UPREGULATES_EXPRESSION, DOWNREGULATES_EXPRESSION,
            new CustomRelationType("acetylates", true),
            new CustomRelationType("deacetylates", true),
            new CustomRelationType("methylates", true),
            new CustomRelationType("demethylates", true),
            new CustomRelationType("activates-gtpase", true),
            new CustomRelationType("inhibits-gtpase", true)
            );

        EdgeSelector ksEdgeSelector = new RelationTypeSelector(PHOSPHORYLATES, DEPHOSPHORYLATES);

        EdgeSelector ppiEdgeSelector = new RelationTypeSelector(INTERACTS_WITH, IN_COMPLEX_WITH);

        String DIR = "/home/ozgunbabur/Analyses/Yo/cis-trans-paths/";

        String[] study = new String[]{"CCRCC", "GBM", "HNSCC", "LSCC", "LUAD", "PDAC", "UCEC"};

        for (int i = 0; i < study.length; i++)
        {
            processInputFile(DIR + "input/" + study[i] + "_tumor_phosphoQTL_gtex_eQTL_coloc_ALL_CIS_GENES.gene_labels.subsetted.tsv",
                cpGraph, cpEdgeSelector, ppiGraph, ppiEdgeSelector, ksEdgeSelector,
                DIR + "networks/" + study[i] + "/", DIR + "output/" + study[i] + "_paths_added.tsv");
        }
    }

    private static void processInputFile(String inFile, SIFGraph cpGraph, EdgeSelector cpEdgeSelector,
                                         SIFGraph ppiGraph, EdgeSelector ppiEdgeSelector, EdgeSelector ksEdgeSelector,
                                         String outDir, String outFile)
    {
        FileUtil.mkdirs(outDir);
        BufferedWriter writer = FileUtil.newBufferedWriter(outFile);

        String[] header = FileUtil.readHeader(inFile);
        int cisInd = ArrayUtil.indexOf(header, "cis_gene_symbol");
        int transInd = ArrayUtil.indexOf(header, "trans_gene_symbol");

        FileUtil.write(header, "\t", writer);
        FileUtil.write("\tOver a kinase path\tCP Distance\tDirected Path\tPPI Distance\tPPI Path", writer);

        int[] max = new int[]{0};
        FileUtil.linesTabbedSkip1(inFile).forEach(t ->
        {
            FileUtil.lnwrite(t, "\t", writer);

            String source = t[cisInd];
            String target = t[transInd];

            Set<Object> kinases = QueryExecutor.searchNeighborhood(cpGraph, ksEdgeSelector, Collections.singleton(target), Direction.UPSTREAM, 1);
            kinases.remove(target);

            Set<Object> interactors = QueryExecutor.searchNeighborhood(ppiGraph, ppiEdgeSelector, Collections.singleton(source), Direction.UNDIRECTED, 1);
            interactors.remove(source);

            Set<Object> common = CollectionUtil.getIntersection(kinases, interactors);

            if (!common.isEmpty())
            {
                FileUtil.tab_write("[" + source + "]-ppi-[" + CollectionUtil.merge(common, ",") + "]-pho-[" + target + "]", writer);
            }
            else
            {
                FileUtil.write("\t", writer);
            }

            Set<Object> result = Collections.emptySet();
            int limit;
            for (limit = 1; result.isEmpty() && limit <= 10; limit++)
            {
                // Run the query
                result = QueryExecutor.searchPathsFromTo(cpGraph, cpEdgeSelector, Collections.singleton(source), Collections.singleton(target), limit);
            }

            if (!result.isEmpty())
            {
                if (limit > max[0]) max[0] = limit;
                System.out.println(source + "\t" + target + "\t" + limit);

                FileUtil.tab_write(limit-1, writer);
                FileUtil.tab_write(toString(result, source, true), writer);

                writeSIF(result, source, target, outDir);
            }
            else
            {
                FileUtil.write("\t\t", writer);
            }

            result = Collections.emptySet();
            for (limit = 1; result.isEmpty() && limit <= 10; limit++)
            {
                // Run the query
                result = QueryExecutor.searchPathsBetween(ppiGraph, ppiEdgeSelector, new HashSet<>(Arrays.asList(source, target)), false, limit);
            }

            if (!result.isEmpty())
            {
                FileUtil.tab_write(limit-1, writer);
                FileUtil.tab_write(toString(result, source, false), writer);
            }
            else
            {
                FileUtil.write("\t\t", writer);
            }
        });

        FileUtil.closeWriter(writer);
        System.out.println("max = " + max[0]);
    }

    private static String toString(Set<Object> graph, String source, boolean directed)
    {
        List<Set<String>> groups = new ArrayList<>();

        Set<String> currentLayer = Collections.singleton(source);

        Set<String> visited = new HashSet<>();
        visited.add(source);

        Set<String> nextLayer = getNextLayer(graph, currentLayer, visited, directed);
        while (!nextLayer.isEmpty())
        {
            groups.add(nextLayer);
            visited.addAll(nextLayer);
            currentLayer = nextLayer;
            nextLayer = getNextLayer(graph, currentLayer, visited, directed);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("[").append(source).append("]");
        for (Set<String> group : groups)
        {
            sb.append("--[").append(CollectionUtil.merge(group, ",")).append("]");
        }
        return sb.toString();
    }

    private static Set<String> getNextLayer(Set<Object> graph, Set<String> currentLayer, Set<String> visited, boolean directed)
    {
        Set<String> nextLayer = new HashSet<>();

        for (Object o : graph)
        {
            if (o instanceof SIFEdge)
            {
                SIFEdge edge = (SIFEdge) o;
                if (currentLayer.contains(edge.getSource()) && !visited.contains(edge.getTarget()))
                {
                    nextLayer.add(edge.getTarget());
                }
                else if (!directed && currentLayer.contains(edge.getTarget()) && !visited.contains(edge.getSource()))
                {
                    nextLayer.add(edge.getSource());
                }
            }
        }

        return nextLayer;
    }

    private static void writeSIF(Set<Object> result, String source, String target, String dir)
    {
        String name = source + "-to-" + target;
        BufferedWriter sifWriter = FileUtil.newBufferedWriter(dir + name + ".sif");

        result.stream().filter(o -> o instanceof SIFEdge).map(o -> (SIFEdge) o).forEach(edge ->
        {
            FileUtil.lnwrite(edge.getSource() + "\t" + edge.getType().getName() + "\t" + edge.getTarget(), sifWriter);
            String mediators = edge.getAnnotation(MEDIATORS);
            if (mediators == null) mediators = "";
            FileUtil.tab_write(mediators, sifWriter);
            String sites = edge.getAnnotation(SITES);
            if (sites != null) FileUtil.tab_write(sites, sifWriter);
        });

        FileUtil.closeWriter(sifWriter);

        BufferedWriter fmtWriter = FileUtil.newBufferedWriter(dir + name + ".format");

        FileUtil.writeln("node\t" + source + "\tbordercolor\t50 50 230", fmtWriter);
        FileUtil.writeln("node\t" + target + "\tbordercolor\t50 50 230", fmtWriter);

        FileUtil.closeWriter(fmtWriter);
    }
}
