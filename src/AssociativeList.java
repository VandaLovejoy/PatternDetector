import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

public class AssociativeList {
    private ArrayList<String []>associativeList;
    private int numbZeros;
    // Constructor
    public AssociativeList() { associativeList = new ArrayList<>();}

    public HashMap<String[], ArrayList> addSpecies(String pathFile, int blockAln, String [] mafTab,
                                                   String RSCAPEBINARY) throws IOException, InterruptedException {

        HashMap<String[], ArrayList> motifs = new HashMap<>();
        String lastDigit = "";
        String finalName = "";
        lastDigit += String.valueOf(blockAln);
        numbZeros = 4 - lastDigit.length();
        int i = 0;
        while (i < numbZeros) {
            finalName += String.valueOf(0);
            i++;
        }
        finalName += lastDigit;

        File file = new File(pathFile +"/stockholm" + "/alifold_" + finalName + ".stk");
        rScape(file, RSCAPEBINARY, pathFile);

        String endFileName = ".fold.sto";
        File dir = new File(pathFile + "/R-Scape");
        File[] matchingFiles = dir.listFiles((dir1, name) -> name.endsWith(endFileName));


        try {
            int v = 0;
            String fileName;
            while (v < matchingFiles.length) {
                fileName = matchingFiles[v].getName();
                String[] arrayName = fileName.split("[_.]");
                int lengthAln = Integer.parseInt(arrayName[4]) - Integer.parseInt(arrayName[3]);
                if (!(lengthAln < 50)) {
                    associativeList = new ArrayList<>();
                    BufferedReader reader = new BufferedReader(new FileReader(matchingFiles[v]));
                    String currentLine = reader.readLine();



                    readBlocks:
                    while (currentLine != null) {
                        if (currentLine.startsWith("#") || currentLine.equals("") || currentLine.startsWith("//")) {
                            currentLine = reader.readLine();
                            continue readBlocks;
                        } else {

                            String[] species = currentLine.split(" ", 2);
                            species[1] = species[1].trim();

                            associativeList.add(species);
                        }

                        currentLine = reader.readLine();
                        continue readBlocks;


                    }

                    int[] cordMotif = getRealCoordinates(Integer.parseInt(arrayName[3]), mafTab);
                    String loci = Arrays.toString(cordMotif);
                    String chrom = mafTab[1].substring((mafTab[1].lastIndexOf(".")) + 1);
                    String lociChrm = chrom + ", " + loci.substring(1, loci.length() - 1) + ", " +
                            mafTab[4] + ", " + arrayName[3] + ", " + arrayName[4];
                    String[] arrayLociChrm = lociChrm.split(", ");
                    motifs.put(arrayLociChrm, associativeList);
                    reader.close();
                }

                v += 1;
            }
             file.delete();



        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        return motifs;
    }

    public int[] getRealCoordinates (int start, String[] mafCord ){
        int nuc = 0;
        int [] cordFinal;
        for (int i = 0; i < start; i++ ){
            char gap = mafCord[6].charAt(i);
            if (!(gap == '-')){
                nuc+=1;
            }
        }
        if (mafCord[4].equals("-")){
            int lociEnd = (Integer.parseInt(mafCord[5] ) + 1  - (Integer.parseInt(mafCord[2]) + nuc)) + 1 ;
            String motifHuman = associativeList.get(0)[1];
            int lociStart = lociEnd - motifHuman.replaceAll("[^ATCGUatcgu]", "").length();
            cordFinal = new int[]{lociStart, lociEnd};
        } else {
            int lociStart = Integer.parseInt(mafCord[2]) + nuc;
            int lociEnd = lociStart + (associativeList.get(0)[1]).length();
            cordFinal = new int[]{lociStart, lociEnd};
        }

           return cordFinal;
    }
    public char[] getAlnTab(String[] seq){
        char[] AlnTab = new char[seq[0].length()];
        if (seq.length>1 && seq[1].length()> 10) {
            AlnTab = seq[1].toCharArray();
        }
        return AlnTab;

    }

    private void rScape(File stockholmFile, String RSCAPEBINARY, String path)throws IOException,
            InterruptedException {
        String pathway = path + "/R-Scape";
        if (!(new File(pathway)).isDirectory())
            (new File(pathway)).mkdirs();

                String cmd = RSCAPEBINARY + " --fold" + " --nofigures "
                + stockholmFile;
        Process process = Runtime.getRuntime().exec(cmd, null, new File(pathway));
        BufferedReader reader =
                new BufferedReader(new InputStreamReader(process.getInputStream()));
        while ((reader.readLine()) != null) {
        }
        process.waitFor();
    }
}
