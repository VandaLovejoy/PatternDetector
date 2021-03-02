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
            dirProgram = "",
            SSZBINARY = "/usr/bin/SISSIz",
            ALIFOLDBINARY = "~/usr/local/bin/RNALalifold",
            RSCAPEBINARY = "/usr/bin/R-scape";

    static double SSZR = -3.0;


    public static synchronized void main(String[] Args) throws IOException, InterruptedException {
        long startTime = System.nanoTime();
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

        GetBinary = Runtime.getRuntime().exec("which SISSIz") ;
        ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));
        if ( (SSZBINARY = ReadBin.readLine() ) == null ) {
            System.out.println("Please install SISSIz (version 2.0), and link it to your $PATH" );
            System.exit(0);
        }

        GetBinary = Runtime.getRuntime().exec("which R-scape") ;
        ReadBin = new BufferedReader(new InputStreamReader(GetBinary.getInputStream()));
        if ( (RSCAPEBINARY = ReadBin.readLine() ) == null ) {
            System.out.println("Please install R-scape, and link it to your $PATH" );
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
                dirProgram = System.getProperty("user.dir");
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

                String cmd = ALIFOLDBINARY + " --id-prefix=alifold" + " --noLP" + " --maxBPspan=200" + " --ribosum_scoring"
                        + " --aln-stk " + Args[Args.length - 1];

               //  executeCommand(cmd);


                ReadFile = new BufferedReader(new FileReader(Args[i]));
                TempTab = new String[lineCount]; // 7 maf columns
                // fill in array from file
                int newLineCount = 0;
                blockAln = 1;
                /************************************************************************
                 ****   This stuff is messy, but avoids problems at last block       ****
                 ************************************************************************/

                ExecutorService MultiThreads = Executors.newFixedThreadPool(NTHREDS);
                List<Future<Runnable>> futures = new ArrayList<Future<Runnable>>();
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

                                ArrayList<String[]> associativeList = new ArrayList<>();
                                mafTabTemp = TempTab[0].split("\\s+");
                                futures = new ArrayList<Future<Runnable>>();
                                // add Path Flag ? < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
                                String Path = OUT_PATH + "/aln/" + mafTabTemp[1].substring(mafTabTemp[1].lastIndexOf(".") + 1);
                                if (!(new File(Path)).isDirectory())
                                    (new File(Path)).mkdirs();
                                int numbZeros = 0;
                                String gcReference = "";
                                String gcSScons = "";

                                HashMap<String[], ArrayList> motifs = new HashMap<>();
                                String lastDigit = "";
                                String finalName = "";
                                lastDigit += String.valueOf(blockAln);
                                numbZeros = 4 - lastDigit.length();
                                int x = 0;
                                while (x < numbZeros) {
                                    finalName += String.valueOf(0);
                                    x++;
                                }
                                finalName += lastDigit;

                                File file = new File(dirProgram + "/" + OUT_PATH + "/stockholm" + "/alifold_"
                                        + finalName + ".stk");


                                try {


                                    BufferedReader reader = new BufferedReader(new FileReader(file));
                                    String currentLine = reader.readLine();

                                    String[] arrayName = new String[5];
                                    int lengthAln = 0;

                                    readAlns:

                                    while (currentLine != null) {


                                        if (currentLine.startsWith("#=GF ID ")) {
                                            arrayName = currentLine.split("[_.]");
                                            lengthAln = Integer.parseInt(arrayName[4]) -
                                                    Integer.parseInt(arrayName[3]);
                                            associativeList = new ArrayList<>();
                                        } else if (currentLine.startsWith("#=GC RF")) {
                                            String[] lineReference = currentLine.split(" ");
                                            gcReference = lineReference[lineReference.length - 1];
                                        } else if (currentLine.startsWith("#=GC SS_cons")) {
                                            String[] lineReference = currentLine.split(" ");
                                            gcSScons = lineReference[lineReference.length - 1];
                                        } else if (currentLine.startsWith("#") || currentLine.equals("")) {
                                            currentLine = reader.readLine();
                                            continue readAlns;

                                        } else if (lengthAln > 50 && !(currentLine.startsWith("//"))) {

                                            String[] species = currentLine.split(" ", 2);
                                            species[1] = species[1].trim();

                                            associativeList.add(species);
                                        }
                                        if ((!associativeList.isEmpty()) && currentLine.startsWith("//")) {
                                            int[] cordMotif = getRealCoordinates(Integer.parseInt(arrayName[3])
                                                    , mafTabTemp, associativeList.get(0)[1]);
                                            String loci = Arrays.toString(cordMotif);
                                            String chrom = mafTabTemp[1].substring((mafTabTemp[1].lastIndexOf(".")) + 1);
                                            String lociChrm = chrom + ", " + loci.substring(1, loci.length() - 1) + ", " +
                                                    mafTabTemp[4] + ", " + arrayName[3] + ", " + arrayName[4] + ", " + gcReference + ", "
                                                    + gcSScons;
                                            String[] arrayLociChrm = lociChrm.split(", ");


                                            Iterator iter = associativeList.iterator();
                                            ArrayList<char[]> alnTab = new ArrayList<>();
                                            while (iter.hasNext()) {
                                                String[] line = (String[]) iter.next();
                                                char[] aln = line[1].toCharArray();
                                                alnTab.add(aln);
                                            }

                                            ScanItFast aln = new ScanItFast( associativeList, alnTab,
                                                    arrayLociChrm,Path,dirProgram + "/" + OUT_PATH, GAPS,
                                                    SSZBINARY, RSCAPEBINARY, VERBOSE, PRINTALL);
                                            aln.setSszR(SSZR);

                                            Future f = MultiThreads.submit(aln);
                                            futures.add(f);
                                        }
                                        currentLine = reader.readLine();
                                        continue readAlns;
                                    }


                                    reader.close();


                                    //  file.delete();


                                } catch (FileNotFoundException e) {
                                    e.printStackTrace();
                                }


                                /*File index = new File(dirProgram + "/"+ OUT_PATH +"/R-Scape");
                                String[] entriesRscape = index.list();
                                for (String s: entriesRscape){
                                    File currentFile = new File(index.getPath(), s);
                                    currentFile.delete();
                                }
*/
                                for (Future<Runnable> g : futures) {
                                    try{
                                    g.get();


                            } catch (Exception Fuck) {
                                System.err.println("MultiThreads took too long!  OOps!");
                                Fuck.printStackTrace();
                            }
                                }


                            blockAln++;
                        } else if (Temp.split("@").length == 2) {
                            blockAln++;
                        }


                        Temp = "";
                    }
                }
            }


                ReadFile.close();
                MultiThreads.shutdown();
                MultiThreads.awaitTermination(60 * 10L, TimeUnit.SECONDS);

            }
        }


        long endTime = System.nanoTime();
        System.out.println("That took"+ (endTime-startTime) +" for program to finish in nanoseconds");
    }

    //*********************************************************************
    //					Alignment generator				*
    //*********************************************************************



    private static int[] getRealCoordinates (int start, String[] mafCord, String motifHuman){
        int [] cordFinal;
        String withoutGap= mafCord[6].substring(0, start);
        int nuc =withoutGap.replaceAll("-", "").length();

        if (mafCord[4].equals("-")){
            int lociEnd = (Integer.parseInt(mafCord[5] ) + 1  - (Integer.parseInt(mafCord[2]) + nuc)) + 1 ;

            int lociStart = lociEnd - motifHuman.replaceAll("[^ATCGUatcgu]", "").length();
            cordFinal = new int[]{lociStart + 1, lociEnd};
        } else {
            int lociStart = Integer.parseInt(mafCord[2]) + nuc;
            int lociEnd = lociStart + (motifHuman.length());
            cordFinal = new int[]{lociStart + 1 , lociEnd};
        }

        return cordFinal;
    }
    protected static char[] getAlnTab(String[] seq){
        char[] AlnTab = new char[seq[0].length()];
        if (seq.length>1 && seq[1].length()> 10) {
            AlnTab = seq[1].toCharArray();
        }
        return AlnTab;

    }



    private static void executeCommand(final String command) throws IOException,
            InterruptedException {

        String Path = dirProgram + "/" + OUT_PATH + "/stockholm";
        if (!(new File(Path)).isDirectory())
            (new File(Path)).mkdirs();
        System.out.println("Executing command " + command);
        Process process = Runtime.getRuntime().exec(command,null, new File(Path));

        BufferedReader reader =
                new BufferedReader(new InputStreamReader(process.getInputStream()));
        while ((reader.readLine()) != null) {
        }
        process.waitFor();

    }

}

