/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package org.lifstools.jmzqc.usecase;

import com.google.common.collect.Range;
import io.github.msdk.MSDKException;
import io.github.msdk.datamodel.ChromatogramType;
import io.github.msdk.io.mzml.MzMLFileImportMethod;
import io.github.msdk.io.mzml.data.MzMLRawDataFile;
import java.io.File;
import java.net.URI;
import java.net.URISyntaxException;
import java.time.OffsetDateTime;
import java.util.AbstractMap.SimpleEntry;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Optional;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.lifstools.jmzqc.AnalysisSoftware;
import org.lifstools.jmzqc.BaseQuality;
import org.lifstools.jmzqc.ControlledVocabulary;
import org.lifstools.jmzqc.CvParameter;
import org.lifstools.jmzqc.InputFile;
import org.lifstools.jmzqc.Metadata;
import org.lifstools.jmzqc.MzQC;
import org.lifstools.jmzqc.QualityMetric;
import org.lifstools.jmzqc.Unit;

/**
 *
 * @author Nils Hoffmann
 */
public class ProteomicsDDAMs1QC {

    private final File inputMzML;

    public ProteomicsDDAMs1QC(File inputMzML) {
        this.inputMzML = inputMzML;
    }

    public static class TicTable {
        List<Float> tic = new ArrayList<>();
        List<Float> rt = new ArrayList<>();
        List<String> nativeSpectrumIdentifier = new ArrayList<>();
        List<Integer> nPeaks = new ArrayList<>();
    }

    public Optional<MzQC> process() throws URISyntaxException {
        MzMLRawDataFile mzMLFile;
        try {
            mzMLFile = new MzMLFileImportMethod(inputMzML.toPath()).execute();
        } catch (MSDKException ex) {
            Logger.getLogger(ProteomicsDDAMs1QC.class.getName()).log(Level.SEVERE, null, ex);
            return Optional.empty();
        }
        var mzMLFormatParameter = new CvParameter("MS:1000584", null, "mzML format", null);
        var inputFile = new InputFile(mzMLFormatParameter, Collections.emptyList(), inputMzML.toURI(), mzMLFile.getName());
        System.out.println("Processing file: " + mzMLFile.getName());
        System.out.println("MS1 mz range...");
        var ms1MzRange = mzMLFile.getScans().stream().filter(
                scan -> scan.getMsLevel() == 1
        ).map(
                scan -> scan.getMzRange()
        ).reduce(
                (l, r) -> l.span(r)
        ).orElse(Range.singleton(Double.NaN));

        System.out.println("TIC and RT values...");

        final TicTable ticTable = new TicTable();
        mzMLFile.getScans().stream().filter(
                scan -> scan.getMsLevel() == 1
        ).forEach(
            scan -> {
                ticTable.tic.add(scan.getTIC());
                ticTable.rt.add(scan.getRetentionTime());
                ticTable.nativeSpectrumIdentifier.add("scan=" + scan.getScanNumber());
                ticTable.nPeaks.add(scan.getNumberOfDataPoints());
            }
        );

        var ms1MzRangeMetric = new QualityMetric("MS:4000069", null, "m/z acquisition range", Arrays.asList(ms1MzRange.lowerEndpoint(), ms1MzRange.upperEndpoint()), null);
        var ticTableMap = new LinkedHashMap<String, List<?>>();
        ticTableMap.put("MS:4000104", ticTable.tic);
        ticTableMap.put("MS:1000894", ticTable.rt);
        ticTableMap.put("MS:1000767", ticTable.nativeSpectrumIdentifier);
        ticTableMap.put("MS:1003059", ticTable.nPeaks);

        var totalIonChromatogram = new QualityMetric(
                "MS:4000104",
                null,
                "total ion currents",
                ticTableMap, null);
        var numberOfChromatogramsMetric = new QualityMetric("MS:4000071", null, "number of chromatograms", mzMLFile.getChromatograms().stream().count(), null);
        System.out.println("RT range...");
        var rtRange = mzMLFile.getScans().stream().map(
                scan -> Range.singleton(scan.getRetentionTime())
        ).reduce(
                (lrt, rrt) -> lrt.span(rrt)
        ).get();

        var rtRangeMetric = new QualityMetric("MS:4000070", null, "retention time acquisition range", Arrays.asList(rtRange.lowerEndpoint(), rtRange.upperEndpoint()), new Unit(new CvParameter("UO:0000010", null, "second", null), null));

        var qualityMetrics = Arrays.asList(
                numberOfChromatogramsMetric,
                ms1MzRangeMetric,
                rtRangeMetric,
                totalIonChromatogram
        );

        var analysisSoftware = new AnalysisSoftware("MS:1000799", null, "custom unreleased software tool", "jmzqc", new URI("https://github.com/MS-Quality-hub/jmzqc"), "1.0.0-RC1");
        Metadata metadata = new Metadata(
                Arrays.asList(analysisSoftware),
                Collections.emptyList(),
                Arrays.asList(inputFile),
                null
        );
        List<BaseQuality> bqs = Arrays.asList(new BaseQuality(metadata, qualityMetrics));
        MzQC mzQC = new MzQC(
                "nils.hoffmann@cebitec.uni-bielefeld.de",
                "Nils Hoffmann",
                Arrays.asList(
                        new ControlledVocabulary(
                                "Proteomics Standards Initiative Mass Spectrometry Ontology",
                                new URI("https://github.com/HUPO-PSI/psi-ms-CV/releases/download/v4.1.103/psi-ms.obo"),
                                "4.1.103"
                        ),
                        new ControlledVocabulary(
                                "Unit Ontology",
                                new URI("https://raw.githubusercontent.com/bio-ontology-research-group/unit-ontology/master/unit.obo"),
                                "f9ff25b"
                        )
                ),
                OffsetDateTime.now(),
                "MzQC for basic TIC QC information",
                bqs,
                Collections.emptyList(),
                "1.0.0");
        return Optional.of(mzQC);
    }

}
