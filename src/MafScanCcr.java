//
//  mafScan.java
//  - - reads a set of maf files, realigns them with mafft, breaks them into windows, calculates stats, 
//  - - scans with SISSIz, outputs bed coordonates of high-confidence predictions 
//  Created by Martin Smith on 27/04/11.
//  Copyright 2013 All rights reserved.
//
import java.util.*; import java.util.concurrent.* ; import java.io.*; import java.math.*; import java.lang.*;
import java.util.Iterator;
public class MafScanCcr {

    static int
            GAPS = 50,
            NTHREDS = 4;
    static boolean
            FILTER_ID = true, //for removing identical sequences
            VERBOSE = false,
            PRINTALL = false,
            REALIGN = false;
    static String
            FILENAME = "",
            OUT_PATH = "",
            SSZBINARY = "/usr/bin/SISSIz",
            ALIFOLDBINARY = "~/usr/local/bin/RNALalifold";

    static double
            SSZ = -2.7,
            SSZR = -2.2;


    public static synchronized void main(String[] Args) throws IOException, InterruptedException {
        // variables
        String[] mafTab;
        String[] mafTabTemp;
        String[] TempTab;
        String Temp = "";
        char[][] AlnTab;
        int blockAln;
        // usage info
        if (Args.length == 0) {
            System.out.println("\n\t\t\t  version 0.9 \n" +
                    "\t __  __          ______    _____  _____          _   _\n" +
                    "\t|  \\/  |   /\\   |  ____|  / ____|/ ____|   /\\   | \\ | |\n" +
                    "\t| \\  / |  /  \\  | |__    | (___ | |       /  \\  |  \\| |\n" +
                    "\t| |\\/| | / /\\ \\ |  __|    \\___ \\| |      / /\\ \\ | . ` |\n" +
                    "\t| |  | |/ ____ \\| |       ____) | |____ / ____ \\| |\\  |\n" +
                    "\t|_|  |_/_/    \\_\\_|      |_____/ \\_____/_/    \\_\\_| \\_|\n" +
                    "\t SCAN MULTIPLE ALIGNMENTS FOR CONSERVED RNA STRUCTURES\n\n" +
                    "Reads a set of maf files, realigns them with mafft, breaks them into windows,\n" +
                    "calculates stats, scans with SISSIz or RNAz, outputs bed coordonates of high-confidence predictions\n\n" +
                    "*** Known issues ***\n" +
                    "-Using MAFFT realignment causes certain jobs to hang when finished;\n this is likely an executor service problem\n\n" +
                    "Usage:     java -jar MafScanCcr.jar [options] -o output/directory -i input.maf (last parameter must be -i)\n\n" +
                    "Output: 	Two types of results are produced:" +
                    "           (1)  the multiple sequence alignments associated to significant predictions \n" +
                    "                are saved to files in the folder specified with the \"-o\" option.\n" +
                    "                File names correspond to their genomic coordinates in a .bed-compatible format. Ex:\n" +
                    "                     output/directory/chrX_12345_12500_80:75:23:14:8:z_300_+.aln\n" +
                    "           (2)  The genomic coordinates (.bed format) of ECSs are also written to the SDOUT\n" +
                    "                (see additional options below)." +
                    "           ***  N.B. the score field corresponds to the SISSIz Z-score x-100, and the RNAz score x100 ***\n\n" +
                    "Options:\n" +
                    "  -all             write out bed entires for all sampled alignments to STDOUT\n" +
                    "  -bs    int       Block Size for splitting large MAF blocks (default 5000)\n" +
                    "  -c     int       number of CPUs for calculations (default 4)\n" +
                    //"  -f               Do not remove identical sequences (not recommended)\n"+
                    "  -g     int       max gap percentage of sequences for 2D prediction (default 50)\n" + // for one or all ??
                    "  -mafft           Realign with mafft-ginsi (slower and buggy, more accurate predictions)\n" +
                    "  -ml    int       Max Length of MAF block for MAFFT realignment (default 10000)\n" +
                    "  -s     int       step size (default 100)\n" +
                    "  -ssz   double    report SISSIz hits below this Z-score (default -2.7)\n" +
                    "  -sszr  double    report SISSIz+RIBOSUM hits below this Z-score (default -2.7)\n" +
                    "  -rnaz  double    report RNAz hits below this SVM probability (default 0.34)\n" +
                    "  -v               verbose (messy but detailed) output\n" +
                    // add a less verbose option, outputing all bed coordinates and scores?
                    "  -w     int       window size (default 200)\n");
            System.exit(0);
        }
        // get binary paths
        Process GetBinary = Runtime.getRuntime().exec("which RNALalifold");
        BufferedReader ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));
        if ((ALIFOLDBINARY = ReadBin.readLine()) == null) {
            System.out.println("Please install RNALalifold and link it to your $PATH");
            System.exit(0);
        }

        ReadBin.close();

        // parse arguments
        for (int i = 0; i != Args.length; i++) {
            if (Args[i].equals("-c")) {  // Threads
                NTHREDS = Integer.parseInt(Args[i + 1]);
                i++;
            } else if (Args[i].equals("-f")) {  // Threads
                FILTER_ID = false;
                i++;
            } else if (Args[i].equals("-g")) {  // gap content
                GAPS = Integer.parseInt(Args[i + 1]);
                i++;
            } else if (Args[i].equals("-o")) { //output directory
                OUT_PATH = Args[i + 1];
                if (!(new File(OUT_PATH)).isDirectory())
                    (new File(OUT_PATH)).mkdirs();
                if (VERBOSE)
                    System.out.println("writing alignments to directory " + OUT_PATH);
                i++;
            } else if (Args[i].equals("-all")) { // step size
                PRINTALL = true;
                i++;
            } else if (Args[i].equals("-v")) { // verbose output
                VERBOSE = true;
            } else if (Args[i].equals("-mafft")) { // realign
                REALIGN = true;
            } else if (Args[i].equals("-i")) {
                i++;
                FILENAME = Args[i].substring(Args[i].lastIndexOf("/") + 1, Args[i].length() - 4);
                //parse out individual alignment blocks from a multi maf file
                int lineCount = 0;
                BufferedReader ReadFile = new BufferedReader(new FileReader(Args[i]));
                String Line = "";
                while ((Line = ReadFile.readLine()) != null)   // count lines for array
                    if (Line.length() > 1 && Line.charAt(0) != '#')
                        lineCount++;

                if (VERBOSE)
                    System.out.println("Read " + (lineCount - 1) + " sequences from file " + FILENAME);
                ReadFile.close();

                /************************************************************************
                 ****   RNALalifold       ****
                 ************************************************************************/

                String cmd = ALIFOLDBINARY + " --id-prefix=alifold" + " --noLP" + " --maxBPspan=200"+ " --ribosum_scoring"
                        + " --aln-stk /scratch/vanda/chromosomes/test";
                executeCommand(cmd);


                ReadFile = new BufferedReader(new FileReader(Args[i]));
                TempTab = new String[lineCount]; // 7 maf columns
                // fill in array from file
                int newLineCount = 0;
                blockAln = 1;
                /************************************************************************
                 ****   This stuff is messy, but avoids problems at last block       ****
                 ************************************************************************/
                readBlocks:
                while ((Line = ReadFile.readLine()) != null || newLineCount == lineCount) {
                    if (Line == null || Line.length() > 1 && Line.charAt(0) != '#')
                        newLineCount++;
                    if (newLineCount > lineCount || Line.length() >= 1) { // only lines with sequences
                        if (newLineCount <= lineCount && Line.length() > 1 && Line.substring(0, 1).equals("s")) {
                            Temp = Temp + Line + "@";
                        } else if ((Line != null && Line.length() >= 1 && Line.substring(0, 1).equals("a")) || newLineCount > lineCount) {
                            if (newLineCount == 1) {

                                continue readBlocks;
                            }

                            if (Temp.split("@").length >= 3) { // at least 3 sequences
                                TempTab = new String[Temp.split("@").length];
                                TempTab = Temp.split("@");
                                Temp = "";
                                AssociativeList species = new AssociativeList();
                                mafTabTemp = TempTab[0].split("\\s+");



                                HashMap<String[], ArrayList> motifs = species.addSpecies(OUT_PATH + "/stockholm/", blockAln, mafTabTemp);
                                for (HashMap.Entry<String[], ArrayList> entry : motifs.entrySet()) {
                                    Iterator iter = entry.getValue().iterator();
                                    ArrayList<char[]> alnTab = new ArrayList<>();
                                    while (iter.hasNext()) {
                                        String[] line = (String[]) iter.next();
                                        char[] aln = species.getAlnTab(line);

                                        alnTab.add(aln);

                                    }

                                    SplitNfold(mafTabTemp, entry.getValue(), alnTab, entry.getKey());


                                }
                                blockAln ++;
                            } else if (Temp.split("@").length == 2) {
                                blockAln ++;
                            }






                              //


                            Temp = "";
                    }
                        }
                    // Wait until all threads are finished
                }

                ReadFile.close();

            }
        }
    }



    //*********************************************************************
    //						scan and parse windows				*
    //*********************************************************************
    private static synchronized boolean SplitNfold(String[] mafTab, ArrayList value, ArrayList<char[]> alnTab,
                                                   String[] key) throws IOException {
        ExecutorService MultiThreads = Executors.newFixedThreadPool(NTHREDS);
        try {
            // check if dir exists
            // add Path Flag ? < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
            String Path = OUT_PATH + "/aln/" + mafTab[1].substring(mafTab[1].lastIndexOf(".") + 1);
            if (!(new File(Path)).isDirectory())
                (new File(Path)).mkdirs();
            List<Future<Runnable>> futures = new ArrayList<Future<Runnable>>();
            ScanItFast Block = new ScanItFast(value, alnTab, key, Path, GAPS, SSZBINARY, VERBOSE, PRINTALL);
             Block.setSsz(SSZ);
             Block.setSszR(SSZR);
            Future f = MultiThreads.submit(Block);
            futures.add(f);
            for (Future<Runnable> g : futures) {
                g.get();
            }
            MultiThreads.shutdown();
            MultiThreads.awaitTermination(60 * 10L, TimeUnit.SECONDS);
        } catch (Exception Fuck) {
            System.err.println("MultiThreads took too long!  OOps!");
            Fuck.printStackTrace();
        }
        return true;
    }


    private static void executeCommand(final String command) throws IOException,
            InterruptedException {

        String Path = OUT_PATH + "/stockholm/";
        if (!(new File(Path)).isDirectory())
            (new File(Path)).mkdirs();
        System.out.println("Executing command " + command);
        Process process = Runtime.getRuntime().exec(command,null, new File("/home/vanda/Downloads/PatternDetector/TEST/stockholm"));

        BufferedReader reader =
                new BufferedReader(new InputStreamReader(process.getInputStream()));
        while ((reader.readLine()) != null) {
        }
        process.waitFor();

    }

}

