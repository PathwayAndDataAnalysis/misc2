package org.panda.misc2;

import org.jetbrains.annotations.NotNull;
import org.panda.resource.siteeffect.Feature;
import org.panda.resource.siteeffect.SiteEffectCollective;
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
//        String klFile = "/home/ozgunbabur/Data/KinaseLibrary/kinase-library-p90-r5.sif";
//        SIFGraph klGraph = cpLoader.load(new FileInputStream(klFile));

//        SiteSpecificGraph klSSGraph = new SiteSpecificGraph("Kinase Library", "phosphorylates");
//        klSSGraph.load(klFile, Collections.singleton("phosphorylates"));

//        klGraph.getAllEdges().forEach(cpGraph::add);

        SiteSpecificGraph phosphorylatesSSGraph = new SiteSpecificGraph("Phospho", "phosphorylates");
//        phosphorylatesSSGraph.load(klFile, Collections.singleton("phosphorylates"));
        phosphorylatesSSGraph.load(cpFile, Collections.singleton("phosphorylates"));
        SiteSpecificGraph dephosphorylatesSSGraph = new SiteSpecificGraph("Dephospho", "dephosphorylates");
        dephosphorylatesSSGraph.load(cpFile, Collections.singleton("dephosphorylates"));
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

        String inputDirName = "/input/";
        String DIR = "/home/ozgunbabur/Analyses/Yo/cis-trans-paths-round-8" + inputDirName;

        Map<String, Map<CausalRelationType, Map<String, Set<String>>>> priorsInMaps = loadCPPriors(cpFile, null);//klFile);

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

                    if (hasCisSlope(filename))
                    {
                        findCausalPaths(filename, outFile.substring(0, outFile.lastIndexOf(".")) + "_CausalPath.tsv", networkDir.substring(0, networkDir.length()-1) + "_CausalPath" + File.separator, priorsInMaps, deduceAnalysisType(filename));
                    }
                }
            }
        });
    }

    private static void writeDirectedSIF(Set<Object> result, String source, String target, double sourceSlope, double targetSlope, String tissue, AnalysisType type, String dir)
    {
        FileUtil.mkdirs(dir);
        String name = source + "-to-" + target;
        if (tissue != null) name += "-in-" + tissue;
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

        FileUtil.writeln("node\t" + source + "\trppasite\t" + source + "_rna|r|" + SLOPE_COLOR.getColorInString(sourceSlope) + "|50 50 50|" + sourceSlope, fmtWriter);

        if (type == AnalysisType.pQTL)
            FileUtil.writeln("node\t" + target + "\tcolor\t" + SLOPE_COLOR.getColorInString(targetSlope), fmtWriter);
        else if (type == AnalysisType.phosphoQTL)
            FileUtil.writeln("node\t" + target + "\trppasite\t" + target + "_site|p|" + SLOPE_COLOR.getColorInString(targetSlope) + "|50 50 50|" + targetSlope, fmtWriter);
        else if (type == AnalysisType.acetylQTL)
            FileUtil.writeln("node\t" + target + "\trppasite\t" + target + "_site|a|" + SLOPE_COLOR.getColorInString(targetSlope) + "|50 50 50|" + targetSlope, fmtWriter);

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
        AnalysisType type = deduceAnalysisType(inFile);

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
            String tissue = cisTissueInd < 0 && transTissueInd < 0 ? null :
                ((cisTissueInd < 0 ? "" : (t[cisTissueInd] + ".")) + t[transTissueInd]).replaceAll(" ", "_");
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
                FileUtil.tab_write("[" + source + "]-cp-[" + CollectionUtil.merge(attachSites(source, common2, phosphoSSGraph, dephosphoSSGraph, acetylSSGraph, deacetylSSGraph), ",") + "]-" + specText + "-[" + target + "]", writer);
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
                FileUtil.tab_write("[" + source + "]-cp-[" + CollectionUtil.merge(attachSites(source, common4, phosphoSSGraph, dephosphoSSGraph, acetylSSGraph, deacetylSSGraph), ",") + "]-" + specText + "-[" + target + "]", writer);
            }
            else FileUtil.write("\t", writer);


            if (!resGraph.isEmpty()) writeDirectedSIF(resGraph, source, target, cisSlope, transSlope, tissue, type, networkDir);
        });

        FileUtil.closeWriter(writer);

        System.out.println("cnt " + Arrays.toString(cnt));

    }

    @NotNull
    private static AnalysisType deduceAnalysisType(String inFile)
    {
        AnalysisType type =
            inFile.contains("phospho") || inFile.contains("kinase") || inFile.contains("phosphatase")  ? AnalysisType.phosphoQTL :
            inFile.contains("acetyl") ? AnalysisType.acetylQTL :
            inFile.contains("pQTL") ? AnalysisType.pQTL : null;
        if (type == null) throw new IllegalArgumentException("Filename not valid: " + inFile);
        return type;
    }

    private static boolean hasCisSlope(String file)
    {
        String[] header = FileUtil.readHeader(file);
        int index = ArrayUtil.indexOf(header, "beta_cis");
        if (index < 0) return false;

        return !FileUtil.linesTabbedSkip1(file).findFirst().get()[index].equals("0");
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

    private static Set<String> attachSites(String source, Set<Object> targets, SiteSpecificGraph... graphs)
    {
        Set<String> withSites = new HashSet<>();
        for (Object target : targets)
        {
            Set<String> sites = new HashSet<>();
            for (SiteSpecificGraph graph : graphs)
            {
                sites.addAll(graph.getSites(source, target.toString()));
            }
            if (sites.isEmpty()) withSites.add(target.toString());
            else
            {
                String name = target + "-" + CollectionUtil.merge(sites, "-");
                withSites.add(name);
            }
        }
        return withSites;
    }

    private static Map<String, Map<CausalRelationType, Map<String, Set<String>>>> loadCPPriors(String cpFile, String klFile)
    {
        Map<String, Map<CausalRelationType, Map<String, Set<String>>>> priors = new HashMap<>();
        fillInCPMap(cpFile, priors);
        if (klFile != null) fillInCPMap(klFile, priors);
        return priors;
    }

    private static void fillInCPMap(String cpFile, Map<String, Map<CausalRelationType, Map<String, Set<String>>>> priors)
    {
        FileUtil.linesTabbed(cpFile).forEach(t ->
        {
            if (!(t[1].endsWith("expression") || t[1].contains("phospho") || t[1].contains("acetyl"))) return;

            String source = t[0];
            if (!priors.containsKey(source)) priors.put(source, new HashMap<>());
            Map<CausalRelationType, Map<String, Set<String>>> relMap = priors.get(source);
            CausalRelationType rel = CausalRelationType.get(t[1]);
            if (!relMap.containsKey(rel)) relMap.put(rel, new HashMap<>());
            Map<String, Set<String>> tarMap = relMap.get(rel);
            String target = t[2];
            if (!tarMap.containsKey(target)) tarMap.put(target, Collections.emptySet());

            if (t.length > 4)
            {
                Set<String> newSites = new HashSet<>(Arrays.asList(t[4].split(";")));
                Set<String> sites = tarMap.get(target);
                if (sites.isEmpty()) tarMap.put(target, newSites);
                else sites.addAll(newSites);
            }
        });
    }

    private static void findCausalPaths(String inputFile, String outFile, String networkDir, Map<String, Map<CausalRelationType, Map<String, Set<String>>>> priors, AnalysisType type)
    {
        String[] header = FileUtil.readHeader(inputFile);
        BufferedWriter writer = FileUtil.newBufferedWriter(outFile);
        SiteEffectCollective sec = new SiteEffectCollective();

        int cisInd = ArrayUtil.indexOf(header, "cis_gene");
        int transInd = ArrayUtil.indexOf(header, "trans_gene");
        int cSlopInd = ArrayUtil.indexOf(header, "beta_cis");
        int tSlopInd = ArrayUtil.indexOf(header, "beta_trans");
        int cisTissueInd = ArrayUtil.indexOf(header, "cis_tissue");
        int transTissueInd = ArrayUtil.indexOf(header, "trans_tissue", "tissue");
        int sitesInd = ArrayUtil.indexOf(header, "PTM_site");

        FileUtil.write(header, "\t", writer);

        FileUtil.linesTabbedSkip1(inputFile).forEach(t ->
        {
            FileUtil.lnwrite(t, "\t", writer);

            String source = t[cisInd];
            String target = t[transInd];
            String tissue = ((cisTissueInd < 0 ? "" : (t[cisTissueInd] + ".")) + t[transTissueInd]).replaceAll(" ", "_");
            double cisSlope = cSlopInd < 0 ? 0 : Double.parseDouble(t[cSlopInd]);
            double transSlope = tSlopInd < 0 ? 0 : Double.parseDouble(t[tSlopInd]);
            Set<String> sites = sitesInd < 0 ? Collections.emptySet() : new HashSet<>(Arrays.asList(t[sitesInd].split("\\|")));

            int dataSign = (int) Math.signum(cisSlope * transSlope);

            Set<String> edges = new HashSet<>();
            List<String> formatLines = new ArrayList<>();
            formatLines.add("node\tall-nodes\tcolor\t255 255 255");
            formatLines.add("node\tall-nodes\tbordercolor\t0 0 0");
            formatLines.add("node\t" + source + "\trppasite\t" + source + "_rna|r|" + SLOPE_COLOR.getColorInString(cisSlope) + "|50 50 50|" + cisSlope);

            Map<CausalRelationType, Map<String, Set<String>>> relMap1 = priors.get(source);
            if (relMap1 != null)
            {
                for (CausalRelationType rel1 : relMap1.keySet())
                {
                    Map<String, Set<String>> midMap = relMap1.get(rel1);
                    for (String middle : midMap.keySet())
                    {
                        // Check if one step search is sufficient
                        if (middle.equals(target) && (type == AnalysisType.pQTL || !CollectionUtil.intersectionEmpty(sites, midMap.get(middle))) && dataSign == rel1.sign)
                        {
                            Set<String> usedSites = CollectionUtil.getIntersection(sites, midMap.get(middle));
                            FileUtil.tab_write(source + " " + rel1 + " " + target + (sites.isEmpty() ? "" : "-" + CollectionUtil.merge(usedSites, "-")), writer);
                            edges.add(source + "\t" + rel1 + "\t" + target);
                            for (String site : usedSites)
                            {
                                Integer effect = sec.getEffect(target, site, rel1.relevantFeature);
                                String line = "node\t" + target + "\trppasite\t" + target + "_" + site + "|p|" + SLOPE_COLOR.getColorInString(transSlope) + "|" + (effect == null || effect == 0 ? "50 50 50" : effect == 1 ? "0 180 20" : "180 0 20") + "|" + transSlope;
                                if (!formatLines.contains(line)) formatLines.add(line);
                            }

                            // No need to go beyond. We found it.
                            continue;
                        } // else go on

                        int siteEffect = 1;
                        Set<String> mSites = new HashSet<>();
                        if (rel1.isPhospho || rel1.isAcetyl)
                        {
                            mSites = midMap.get(middle);

                            Integer effect = sec.getEffect(middle, mSites, rel1.relevantFeature);
                            siteEffect = Objects.requireNonNullElse(effect, 0);
                        }

                        if (siteEffect != 0)
                        {
                            Map<CausalRelationType, Map<String, Set<String>>> relMap2 = priors.get(middle);
                            if (relMap2 != null)
                            {
                                for (CausalRelationType rel2 : relMap2.keySet())
                                {
                                    if ((type ==AnalysisType.pQTL && rel2.isExpression) ||
                                        (type ==AnalysisType.phosphoQTL && rel2.isPhospho) ||
                                        (type ==AnalysisType.acetylQTL && rel2.isAcetyl))
                                    {
                                        Map<String, Set<String>> tarMap = relMap2.get(rel2);
                                        if (tarMap.containsKey(target))
                                        {
                                            if (type == AnalysisType.pQTL || !CollectionUtil.intersectionEmpty(tarMap.get(target), sites))
                                            {
                                                if (dataSign == rel1.sign * siteEffect * rel2.sign)
                                                {
                                                    Set<String> usedSites = CollectionUtil.getIntersection(sites, tarMap.get(target));
                                                    FileUtil.tab_write(source + " " + rel1 + " " + middle +
                                                        (mSites.isEmpty() ? "" : "-" + CollectionUtil.merge(mSites, "-") + (siteEffect == 1 ? "(a)" : "(i)")) +
                                                        " " + rel2 + " " + target + (sites.isEmpty() ? "" : "-" + CollectionUtil.merge(usedSites, "-")), writer);

                                                    edges.add(source + "\t" + rel1 + "\t" + middle);
                                                    edges.add(middle + "\t" + rel2 + "\t" + target);
                                                    for (String mSite : mSites)
                                                    {
                                                        Integer effect = sec.getEffect(middle, mSite, rel1.relevantFeature);
                                                        if (effect != null && effect != 0)
                                                            formatLines.add("node\t" + middle + "\trppasite\t" + middle + "_" + mSite + "|p|255 255 255|" + (effect == 1 ? "0 180 20" : "180 0 20") + "|0");
                                                    }
                                                    for (String site : usedSites)
                                                    {
                                                        Integer effect = sec.getEffect(target, site, rel2.relevantFeature);
                                                        String line = "node\t" + target + "\trppasite\t" + target + "_" + site + "|p|" + SLOPE_COLOR.getColorInString(transSlope) + "|" + (effect == null || effect == 0 ? "50 50 50" : effect == 1 ? "0 180 20" : "180 0 20") + "|" + transSlope;
                                                        if (!formatLines.contains(line)) formatLines.add(line);
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            if (!edges.isEmpty())
            {
                FileUtil.createDirectories(networkDir);
                String name = source + "-to-" + target + "-in-" + tissue + "-CausalPath";
                FileUtil.writeLinesToFile(edges, networkDir + name + ".sif");
                FileUtil.writeLinesToFile(formatLines, networkDir + name + ".format");
            }
        });
        FileUtil.closeWriter(writer);
    }
    
    private enum AnalysisType
    {
        pQTL, phosphoQTL, acetylQTL
    }

    private enum CausalRelationType
    {
        phosphorylates(1, Feature.PHOSPHORYLATION, true, false, false),
        dephosphorylates(-1, Feature.PHOSPHORYLATION, true, false, false),
        upregulates_expression(1, Feature.GLOBAL_PROTEIN, false, false, true),
        downregulates_expression(-1, Feature.GLOBAL_PROTEIN, false, false, true),
        acetylates(1, Feature.ACETYLATION, false, true, false),
        deacetylates(-1, Feature.ACETYLATION, false, true, false);

        CausalRelationType(int sign, Feature relevantFeature, boolean isPhospho, boolean isAcetyl, boolean isExpression)
        {
            this.sign = sign;
            this.relevantFeature = relevantFeature;
            this.isPhospho = isPhospho;
            this.isAcetyl = isAcetyl;
            this.isExpression = isExpression;
        }

        int sign;
        Feature relevantFeature;
        boolean isPhospho;
        boolean isAcetyl;
        boolean isExpression;

        @Override
        public String toString()
        {
            return super.toString().replaceAll("_", "-");
        }

        public static CausalRelationType get(String name)
        {
            return CausalRelationType.valueOf(name.replaceAll("-", "_"));
        }
    }
}
