package org.panda.misc2;

import org.panda.utility.ArrayUtil;
import org.panda.utility.CollectionUtil;
import org.panda.utility.FileUtil;
import org.panda.utility.ValToColor;
import org.panda.utility.graph.SiteSpecificGraph;
import org.pathwaycommons.sif.io.Loader;
import org.pathwaycommons.sif.model.CustomRelationType;
import org.pathwaycommons.sif.model.SIFEdge;
import org.pathwaycommons.sif.model.SIFGraph;
import org.pathwaycommons.sif.query.Direction;
import org.pathwaycommons.sif.query.QueryExecutor;
import org.pathwaycommons.sif.util.EdgeAnnotationType;
import org.pathwaycommons.sif.util.EdgeSelector;
import org.pathwaycommons.sif.util.RelationTypeSelector;

import java.awt.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.*;
import java.util.List;
import java.util.stream.Collectors;

import static org.pathwaycommons.sif.model.RelationTypeEnum.INTERACTS_WITH;
import static org.pathwaycommons.sif.model.RelationTypeEnum.IN_COMPLEX_WITH;
import static org.pathwaycommons.sif.model.SignedTypeEnum.*;
import static org.pathwaycommons.sif.util.EdgeAnnotationType.*;

public class PathFinder
{
    static ValToColor SLOPE_COLOR = new ValToColor(new double[]{-1.5, 0, 1.5}, new Color[]{new Color(80, 80, 200), Color.white, new Color(200, 80, 80)});

    public static void main(String[] args) throws IOException {
        // Tell the type and order of edge annotations in the SIF resource
        EdgeAnnotationType[] cpAnnotTypes = new EdgeAnnotationType[]{MEDIATORS, SITES};
        EdgeAnnotationType[] ppiAnnotTypes = new EdgeAnnotationType[]{};

        // Initialize loader
        Loader cpLoader = new Loader(cpAnnotTypes);
        Loader ppiLoader = new Loader(ppiAnnotTypes);

        // Load graphs
        String cpFile = "/home/ozgunbabur/Documents/causal-priors.txt";
        SIFGraph cpGraph = cpLoader.load(new FileInputStream(cpFile));
        String klFile = "/home/ozgunbabur/Data/KinaseLibrary/kinase-library-p90-r15.sif";
        SIFGraph klGraph = cpLoader.load(new FileInputStream(klFile));

        SiteSpecificGraph klSSGraph = new SiteSpecificGraph("Kinase Library", "phosphorylates");
        klSSGraph.load(klFile, Collections.singleton("phosphorylates"));

        klGraph.getAllEdges().forEach(cpGraph::add);

        SiteSpecificGraph phosphorylatesSSGraph = new SiteSpecificGraph("Phospho", "phosphorylates");
        phosphorylatesSSGraph.load(klFile, Collections.singleton("phosphorylates"));
        phosphorylatesSSGraph.load(cpFile, Collections.singleton("phosphorylates"));
        SiteSpecificGraph dephosphorylatesSSGraph = new SiteSpecificGraph("Dephospho", "dephosphorylates");
        phosphorylatesSSGraph.load(cpFile, Collections.singleton("phosphorylates"));
        SiteSpecificGraph acetylatesSSGraph = new SiteSpecificGraph("Acetyl", "acetylates");
        acetylatesSSGraph.load(cpFile, Collections.singleton("acetylates"));
        SiteSpecificGraph deacetylatesSSGraph = new SiteSpecificGraph("Deacetyl", "deacetylates");
        deacetylatesSSGraph.load(cpFile, Collections.singleton("deacetylates"));

        SIFGraph ppiGraph = ppiLoader.load(new FileInputStream("/home/ozgunbabur/Documents/PathwayCommons12.All.hgnc.sif"));

        //--- Perform a query on the graph

        // The query will traverse only two type of relations in this example
        CustomRelationType acetylatesRelationType = new CustomRelationType("acetylates", true);
        CustomRelationType deacetylatesRelationType = new CustomRelationType("deacetylates", true);
        EdgeSelector cpEdgeSelector = new RelationTypeSelector(PHOSPHORYLATES, DEPHOSPHORYLATES, UPREGULATES_EXPRESSION, DOWNREGULATES_EXPRESSION,
            acetylatesRelationType,
            deacetylatesRelationType,
            new CustomRelationType("methylates", true),
            new CustomRelationType("demethylates", true),
            new CustomRelationType("activates-gtpase", true),
            new CustomRelationType("inhibits-gtpase", true)
            );

        EdgeSelector ksEdgeSelector = new RelationTypeSelector(PHOSPHORYLATES, DEPHOSPHORYLATES);
        EdgeSelector tfEdgeSelector = new RelationTypeSelector(UPREGULATES_EXPRESSION, DOWNREGULATES_EXPRESSION);
        EdgeSelector acetylEdgeSelector = new RelationTypeSelector(acetylatesRelationType, deacetylatesRelationType);

        EdgeSelector ppiEdgeSelector = new RelationTypeSelector(INTERACTS_WITH, IN_COMPLEX_WITH);

        String inputDirName = "/rand-input/";
        String DIR = "/home/ozgunbabur/Analyses/Yo/cis-trans-paths-round-4" + inputDirName;

        FileUtil.processDirsRecursive(new File(DIR), dir ->
        {
            for (File file : dir.listFiles())
            {
                String filename = file.getPath();
                if (filename.endsWith(".tsv"))
                {
                    String networkDir = filename.replaceFirst(inputDirName, "/networks/");
                    networkDir = networkDir.substring(0, networkDir.lastIndexOf("/") + 1);
                    String anaName = filename.substring(filename.lastIndexOf(File.separator) + 1, filename.lastIndexOf("."));
                    anaName = anaName.replaceAll("causalpath_input\\.", "");
                    networkDir += anaName + File.separator;


                    String outFile = filename.replaceFirst(inputDirName, "/output/");
//                    processInputFile(filename,
//                        cpGraph, cpEdgeSelector, ppiGraph, ppiEdgeSelector, ksEdgeSelector, acetylEdgeSelector, tfEdgeSelector
//                        networkDir, outFile);
                    processInputFileTFOnly(filename, cpGraph, tfEdgeSelector, ksEdgeSelector, acetylEdgeSelector, cpEdgeSelector, ppiGraph, ppiEdgeSelector, phosphorylatesSSGraph, dephosphorylatesSSGraph, acetylatesSSGraph, deacetylatesSSGraph, networkDir, outFile);
                }
            }
        });
    }

    private static void processInputFile(String inFile, SIFGraph cpGraph, EdgeSelector cpEdgeSelector,
                                         SIFGraph ppiGraph, EdgeSelector ppiEdgeSelector, EdgeSelector ksEdgeSelector,
                                         EdgeSelector acetylEdgeSelector, EdgeSelector tfEdgeSelector, String networkDir, String outFile)
    {
        FileUtil.mkdirsOfFilePath(outFile);

        BufferedWriter writer = FileUtil.newBufferedWriter(outFile);

        boolean acetyl = inFile.contains("acetyl");

        String[] header = FileUtil.readHeader(inFile);
        int cisInd = ArrayUtil.indexOf(header, "cis_gene");
        int transInd = ArrayUtil.indexOf(header, "trans_gene");
        int cSlopInd = ArrayUtil.indexOf(header, "cis_slope");
        int tSlopInd = ArrayUtil.indexOf(header, "trans_slope");
        int tissueInd = ArrayUtil.indexOf(header, "tissue");

        FileUtil.write(header, "\t", writer);
        FileUtil.write("\tOver "+ (acetyl ? "an acetyl transferase" : "a kinase") + " path\tCP Distance\tDirected Path\tPPI Distance\tPPI Path", writer);

        int[] max = new int[]{0};
        FileUtil.linesTabbedSkip1(inFile).forEach(t ->
        {
            FileUtil.lnwrite(t, "\t", writer);

            String source = t[cisInd];
            String target = t[transInd];

            double cisSlope = Double.parseDouble(t[cSlopInd]);
            double transSlope = Double.parseDouble(t[tSlopInd]);

            String tissue = t[tissueInd];

            Set<Object> mids = QueryExecutor.searchNeighborhood(cpGraph, acetyl ? acetylEdgeSelector : ksEdgeSelector, Collections.singleton(target), Direction.UPSTREAM, 1);
            mids.remove(target);

            Set<Object> interactors = QueryExecutor.searchNeighborhood(ppiGraph, ppiEdgeSelector, Collections.singleton(source), Direction.UNDIRECTED, 1);
            interactors.remove(source);

            Set<Object> common = CollectionUtil.getIntersection(mids, interactors);

            if (!common.isEmpty())
            {
                FileUtil.tab_write("[" + source + "]-ppi-[" + CollectionUtil.merge(common, ",") + "]-pho-[" + target + "]", writer);
                writeOverAKinaseSIF(source, target, common.stream().map(Object::toString).collect(Collectors.toSet()), cisSlope, transSlope, tissue, networkDir.replaceFirst("/networks/", "/networks/over-" + (acetyl ? "an-acetyl-transferase" : "a-kinase") + "/"));
            }
            else
            {
                FileUtil.write("\t", writer);
            }

            Set<Object> result = Collections.emptySet();
            int limit;
            for (limit = 1; result.isEmpty() && limit <= 3; limit++)
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

                writeDirectedSIF(result, source, target, cisSlope, transSlope, tissue, networkDir.replaceFirst("/networks/", "/networks/directed/"));
            }
            else
            {
                FileUtil.write("\t\t", writer);
            }

            result = Collections.emptySet();
            for (limit = 1; result.isEmpty() && limit <= 1; limit++)
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

    private static void writeDirectedSIF(Set<Object> result, String source, String target, double sourceSlope, double targetSlope, String tissue, String dir)
    {
        FileUtil.mkdirs(dir);
        String name = source + "-to-" + target + "-in-" + tissue;
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

        FileUtil.writeln("node\t" + source + "\tcolor\t" + SLOPE_COLOR.getColorInString(sourceSlope), fmtWriter);
        FileUtil.writeln("node\t" + target + "\tcolor\t" + SLOPE_COLOR.getColorInString(targetSlope), fmtWriter);

        FileUtil.closeWriter(fmtWriter);
    }

    private static void writeOverAKinaseSIF(String source, String target, Set<String> kinases, double sourceSlope, double targetSlope, String tissue, String dir)
    {
        FileUtil.mkdirs(dir);
        String name = source + "-to-" + target + "-in-" + tissue;
        BufferedWriter sifWriter = FileUtil.newBufferedWriter(dir + name + ".sif");

        kinases.forEach(kinase ->
        {
            FileUtil.lnwrite( source + "\tinteracts-with\t" + kinase, sifWriter);
            FileUtil.lnwrite( kinase + "\tphosphorylates\t" + target, sifWriter);


        });

        FileUtil.closeWriter(sifWriter);

        BufferedWriter fmtWriter = FileUtil.newBufferedWriter(dir + name + ".format");

        FileUtil.writeln("node\t" + source + "\tcolor\t" + SLOPE_COLOR.getColorInString(sourceSlope), fmtWriter);
        FileUtil.writeln("node\t" + target + "\tcolor\t" + SLOPE_COLOR.getColorInString(targetSlope), fmtWriter);

        FileUtil.closeWriter(fmtWriter);
    }

    static void processInputFileTFOnly(String inFile, SIFGraph cpGraph, EdgeSelector tfEdgeSelector,
                                       EdgeSelector ksEdgeSelector, EdgeSelector acetylEdgeSelector,
                                       EdgeSelector cpEdgeSelector, SIFGraph ppiGraph, EdgeSelector ppiEdgeSelector,
                                       SiteSpecificGraph phosphoSSGraph, SiteSpecificGraph dephosphoSSGraph,
                                       SiteSpecificGraph acetylSSGraph, SiteSpecificGraph deacetylSSGraph,
                                       String networkDir, String outFile)
    {
        System.out.println("inFile = " + inFile);
        FileUtil.mkdirsOfFilePath(outFile);
        AnalysisType type = inFile.contains("phosphoQTL") ? AnalysisType.phosphoQTL : inFile.contains("acetylQTL") ? AnalysisType.acetylQTL : inFile.contains("pQTL") ? AnalysisType.pQTL : null;
        if (type == null) throw new IllegalArgumentException("Filename not valid: " + inFile);

        BufferedWriter writer = FileUtil.newBufferedWriter(outFile);

        String[] header = FileUtil.readHeader(inFile);
        int cisInd = ArrayUtil.indexOf(header, "cis_gene");
        int transInd = ArrayUtil.indexOf(header, "trans_gene");
        int cSlopInd = ArrayUtil.indexOf(header, "beta_cis");
        int tSlopInd = ArrayUtil.indexOf(header, "beta_trans");
        int cisTissueInd = ArrayUtil.indexOf(header, "cis_tissue");
        int transTissueInd = ArrayUtil.indexOf(header, "trans_tissue", "tissue");
        int sitesInd = ArrayUtil.indexOf(header, "PTM_site");

        FileUtil.write(header, "\t", writer);
        FileUtil.write(type == AnalysisType.pQTL ? "\tDirect expression link\tDirect CP link\tDirect PPI link\tPPI - TF path\tCP - TF path" :
            type == AnalysisType.phosphoQTL ? "\tDirect phospho link\tDirect phospho link site-specific\tDirect CP link\tDirect PPI link\tPPI - phos path\tCP - phos path\tPPI - phos path site-specific\tCP - phos path site-specific" :
            "\tDirect acetyl link\tDirect acetyl link site-specific\tDirect CP link\tDirect PPI link\tPPI - acetyl path\tCP - acetyl path\tPPI - acetyl path site-specific\tCP - acetyl path site-specific", writer);

        int[] cnt = new int[]{0, 0, 0, 0};
        FileUtil.linesTabbedSkip1(inFile).forEach(t ->
        {
            FileUtil.lnwrite(t, "\t", writer);

            String source = t[cisInd];
            String target = t[transInd];
            String tissue = ((cisTissueInd < 0 ? "" : (t[cisTissueInd] + ".")) + t[transTissueInd]).replaceAll(" ", "_");
            double cisSlope = cSlopInd < 0 ? 0 : Double.parseDouble(t[cSlopInd]);
            double transSlope = tSlopInd < 0 ? 0 : Double.parseDouble(t[tSlopInd]);
            Set<String> sites = sitesInd < 0 ? Collections.emptySet() : new HashSet<>(Arrays.asList(t[sitesInd].split("\\|")));

            Set<Object> resGraph = new HashSet<>();

            Set<Object> specUpstrQueryResult = QueryExecutor.searchNeighborhood(cpGraph, type == AnalysisType.phosphoQTL ? ksEdgeSelector : type == AnalysisType.acetylQTL ? acetylEdgeSelector : tfEdgeSelector, Collections.singleton(target), Direction.UPSTREAM, 1);
            Set<SIFEdge> specUpstreamEdges = specUpstrQueryResult.stream().filter(o -> o instanceof SIFEdge).map(o -> (SIFEdge) o).collect(Collectors.toSet());
            Set<SIFEdge> directSpecEdges = specUpstreamEdges.stream().filter(e -> e.getSource().equals(source)).collect(Collectors.toSet());

            Set<Object> cpDwstrQueryResult = QueryExecutor.searchNeighborhood(cpGraph, cpEdgeSelector, Collections.singleton(source), Direction.DOWNSTREAM, 1);

            Set<String> ssUstrGenes = new HashSet<>();
            if (type == AnalysisType.phosphoQTL)
            {
                ssUstrGenes.addAll(phosphoSSGraph.getUpstream(target, sites));
                ssUstrGenes.addAll(dephosphoSSGraph.getUpstream(target, sites));
            }
            else if (type == AnalysisType.acetylQTL)
            {
                ssUstrGenes.addAll(acetylSSGraph.getUpstream(target, sites));
                ssUstrGenes.addAll(deacetylSSGraph.getUpstream(target, sites));
            }

            // Direct specific edges
            if (directSpecEdges.size() > 1) System.err.println("WARNING: bipolar relations from " + source + " to " + target);
            String priorCell = "";
            if (!directSpecEdges.isEmpty())
            {
                SIFEdge edge = directSpecEdges.iterator().next();
                resGraph.add(edge);
                FileUtil.tab_write(edge.toString(), writer);
                priorCell = edge.toString();
                cnt[0]++;
            }
            else FileUtil.tab_write("", writer);

            // Direct site-spec edges
            if (type != AnalysisType.pQTL)
            {
                FileUtil.tab_write(ssUstrGenes.contains(source) ? priorCell : "", writer);
            }

            // Direct CP edges
            if (cpDwstrQueryResult.contains(target))
            {
                Set<SIFEdge> directCPEdges = cpDwstrQueryResult.stream().filter(o -> o instanceof SIFEdge).map(o -> (SIFEdge) o).filter(e -> e.getSource().equals(source) && e.getTarget().equals(target)).collect(Collectors.toSet());
                if (directCPEdges.isEmpty()) throw new RuntimeException("Cannot be empty!!");
                resGraph.addAll(directCPEdges);
                SIFEdge edge = directCPEdges.iterator().next();
                FileUtil.tab_write(edge.toString(), writer);
            }
            else FileUtil.tab_write("", writer);

            Set<Object> midGenes = specUpstreamEdges.stream().map(SIFEdge::getSource).collect(Collectors.toSet());

            Set<Object> interactors = QueryExecutor.searchNeighborhood(ppiGraph, ppiEdgeSelector, Collections.singleton(source), Direction.UNDIRECTED, 1);
            interactors.remove(source);

            // Direct PPI interactors

            boolean directPPI = interactors.contains(target);
            if (directPPI)
            {
                FileUtil.tab_write(source + "-ppi-" + target, writer);
                putPPIsInResults(resGraph, source, target, interactors);
                cnt[1]++;
            }
            else FileUtil.write("\t", writer);

            // PPI - spec path

            Set<Object> common = CollectionUtil.getIntersection(midGenes, interactors);

            String specText = type == AnalysisType.pQTL ? "exp" : type == AnalysisType.phosphoQTL ? "phos" : "acet";
            if (!common.isEmpty())
            {
                FileUtil.tab_write("[" + source + "]-ppi-[" + CollectionUtil.merge(common, ",") + "]-" + specText + "-[" + target + "]", writer);
                putPPIsInResults(resGraph, source, common, interactors);
                specUpstreamEdges.stream().filter(e -> common.contains(e.getSource())).forEach(resGraph::add);
                cnt[2]++;
            }
            else FileUtil.write("\t", writer);

            // CP - spec path

            Set<Object> common2 = CollectionUtil.getIntersection(midGenes, cpDwstrQueryResult);

            if (!common2.isEmpty())
            {
                FileUtil.tab_write("[" + source + "]-cp-[" + CollectionUtil.merge(common2, ",") + "]-" + specText + "-[" + target + "]", writer);
                cpDwstrQueryResult.stream().filter(o -> o instanceof SIFEdge).map(o -> (SIFEdge) o).filter(e -> common2.contains(e.getTarget())).forEach(resGraph::add);
                specUpstreamEdges.stream().filter(e -> common2.contains(e.getSource())).forEach(resGraph::add);
                cnt[3]++;
            }
            else FileUtil.write("\t", writer);

            // PPI - spec site specific

            Set<Object> common3 = CollectionUtil.getIntersection(new HashSet<>(ssUstrGenes), interactors);
            if (!common3.isEmpty())
            {
                FileUtil.tab_write("[" + source + "]-ppi-[" + CollectionUtil.merge(common3, ",") + "]-" + specText + "-[" + target + "]", writer);
            }
            else FileUtil.write("\t", writer);

            // CP - spec site specific
            Set<Object> common4 = CollectionUtil.getIntersection(new HashSet<>(ssUstrGenes), cpDwstrQueryResult);
            if (!common4.isEmpty())
            {
                FileUtil.tab_write("[" + source + "]-cp-[" + CollectionUtil.merge(common4, ",") + "]-" + specText + "-[" + target + "]", writer);
            }
            else FileUtil.write("\t", writer);


            if (!resGraph.isEmpty()) writeDirectedSIF(resGraph, source, target, cisSlope, transSlope, tissue, networkDir);
        });

        FileUtil.closeWriter(writer);

        System.out.println("cnt " + Arrays.toString(cnt));

    }

    private static void putPPIsInResults(Set<Object> resGraph, String source, String target, Set<Object> interactors)
    {
        interactors.stream().filter(o -> o instanceof SIFEdge).map(o -> (SIFEdge)o)
            .filter(e -> (e.getSource().equals(source) && e.getTarget().equals(target)) || (e.getSource().equals(target) && e.getTarget().equals(source)))
            .forEach(resGraph::add);
    }
    private static void putPPIsInResults(Set<Object> resGraph, String gene1, Set<Object> others, Set<Object> interactors)
    {
        interactors.stream().filter(o -> o instanceof SIFEdge).map(o -> (SIFEdge)o)
            .filter(e -> (e.getSource().equals(gene1) && others.contains(e.getTarget()) || (others.contains(e.getSource()) && e.getTarget().equals(gene1))))
            .forEach(resGraph::add);
    }

    private enum AnalysisType
    {
        pQTL, phosphoQTL, acetylQTL
    }
}
