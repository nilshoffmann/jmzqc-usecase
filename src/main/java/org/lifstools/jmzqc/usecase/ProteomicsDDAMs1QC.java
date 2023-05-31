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

        var tic = mzMLFile.getScans().stream().filter(
                scan -> scan.getMsLevel() == 1
        ).map(
                scan -> scan.getTIC()
        ).toList();

        var rt = mzMLFile.getScans().stream().filter(
                scan -> scan.getMsLevel() == 1
        ).map(
                scan -> scan.getRetentionTime()
        ).toList();

        var nativeSpectrumIdentifier = mzMLFile.getScans().stream().filter(
                scan -> scan.getMsLevel() == 1
        ).map(
                scan -> "scan=" + scan.getScanNumber()
        ).toList();

        var nPeaks = mzMLFile.getScans().stream().filter(
                scan -> scan.getMsLevel() == 1
        ).map(
                scan -> scan.getNumberOfDataPoints()
        ).toList();

        var ms1MzRangeMetric = new QualityMetric("MS:4000069", null, "m/z acquisition range", Arrays.asList(ms1MzRange.lowerEndpoint(), ms1MzRange.upperEndpoint()), null);
        System.out.println("TIC and RT values...");
        var ticValuesAndRts = mzMLFile.getChromatograms().stream().filter(
                chrom -> chrom.getChromatogramType() == ChromatogramType.TIC
        ).findFirst().map(
                chrom -> {
                    return new SimpleEntry<>(chrom.getRetentionTimes(), chrom.getIntensityValues());
                }
        ).orElse(new SimpleEntry<>(new float[0], new float[0]));

        var ticTable = new LinkedHashMap<String, List<?>>();
        ticTable.put("MS:4000104", tic);
        ticTable.put("MS:1000894", rt);
        ticTable.put("MS:1000767", nativeSpectrumIdentifier);
        ticTable.put("MS:1003059", nPeaks);

        var totalIonChromatogram = new QualityMetric(
                "MS:4000104",
                null,
                "total ion currents",
                ticTable, null);
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
