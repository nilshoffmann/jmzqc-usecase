package org.lifstools.jmzqc.usecase;

import com.fasterxml.jackson.annotation.JsonInclude.Include;
import com.fasterxml.jackson.core.JsonFactoryBuilder;
import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.core.json.JsonReadFeature;
import com.fasterxml.jackson.databind.DeserializationContext;
import com.fasterxml.jackson.databind.JsonDeserializer;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectWriter;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.databind.module.SimpleModule;
import com.google.common.collect.Range;
import com.networknt.schema.ValidationMessage;
import io.github.msdk.MSDKException;
import io.github.msdk.MSDKRuntimeException;
import io.github.msdk.datamodel.ChromatogramType;
import io.github.msdk.io.mzml.MzMLFileImportMethod;
import io.github.msdk.io.mzml.data.MzMLRawDataFile;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.StandardCopyOption;
import java.time.OffsetDateTime;
import java.util.AbstractMap.SimpleEntry;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.knowm.xchart.BitmapEncoder;
import org.knowm.xchart.QuickChart;
import org.knowm.xchart.XYChart;
import org.lifstools.jmzqc.AnalysisSoftware;
import org.lifstools.jmzqc.BaseQuality;
import org.lifstools.jmzqc.ControlledVocabulary;
import org.lifstools.jmzqc.Converter;
import org.lifstools.jmzqc.Coordinate;
import org.lifstools.jmzqc.CvParameter;
import org.lifstools.jmzqc.InputFile;
import org.lifstools.jmzqc.Metadata;
import org.lifstools.jmzqc.MzQC;
import org.lifstools.jmzqc.QualityMetric;
import org.lifstools.jmzqc.Unit;

/**
 *
 * @author nilshoffmann
 */
public class JmzqcProteomics {

    public static void main(String[] args) {
        var outputDir = new File("proteomics-qc");
        var baseName = "20181113_010_autoQC01";
        var file = new File(outputDir, baseName + ".mzML");
        var fileUrl = "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000086542/ccms_peak/20181113_010_autoQC01.mzML&forceDownload=true";

        try {
            if (outputDir.isDirectory()) {
                // clear all files if they exist in target dir
                Files.walk(outputDir.toPath())
                        .sorted(Comparator.reverseOrder())
                        .map(Path::toFile)
                        .forEach(File::delete);
            }
            outputDir.mkdirs();
            System.out.println("Downloading file");
            InputStream in = new URL(fileUrl).openStream();
            Files.copy(in, file.toPath(), StandardCopyOption.REPLACE_EXISTING);
            System.out.println("Done!");
        } catch (IOException ex) {
            System.err.println("Exception:" + ex.getLocalizedMessage());
        }

        List<MzMLRawDataFile> mzMLData = Collections.emptyList();
        try {
            var mzMLFilePaths = Files.list(outputDir.toPath()).collect(Collectors.toList());
            mzMLData = mzMLFilePaths.stream().map(path -> {
                try {
                    return new MzMLFileImportMethod(path).execute();
                } catch (MSDKException ex) {
                    throw new MSDKRuntimeException(ex);
                }
            }
            ).collect(Collectors.toList());
        } catch (IOException ex) {
            System.err.println("Exception:" + ex.getLocalizedMessage());
        }
        var mzMLFormatParameter = new CvParameter("MS:1000584", null, "mzML format", null);
        Map<InputFile, List<QualityMetric>> mzMLFileStats = mzMLData.stream().map((t) -> {
            System.out.println("Processing file: " + t.getName());
            System.out.println("MS1 mz range...");
            var ms1MzRange = t.getScans().stream().filter(
                    scan -> scan.getMsLevel() == 1
            ).map(
                    scan -> scan.getMzRange()
            ).reduce(
                    (l, r) -> l.span(r)
            ).orElse(Range.singleton(Double.NaN));

            var ms1MzRangeMetric = new QualityMetric("MS:4000069", null, "m/z acquisition range", Arrays.asList(ms1MzRange.lowerEndpoint(), ms1MzRange.upperEndpoint()), null);
            System.out.println("TIC and RT values...");
            var ticValuesAndRts = t.getChromatograms().stream().filter(
                    chrom -> chrom.getChromatogramType() == ChromatogramType.TIC
            ).findFirst().map(
                    chrom -> {
                        return new SimpleEntry<>(chrom.getRetentionTimes(), chrom.getIntensityValues());
                    }
            ).orElse(new SimpleEntry<>(new float[0], new float[0]));
            var totalIonChromatogram = new QualityMetric(
                    "MS:1000235",
                    null,
                    "total ion current chromatogram",
                    ticValuesAndRts.getValue(), new Unit(new CvParameter("UO:0000010", null, "second", ticValuesAndRts.getKey()), null));
            var numberOfChromatogramsMetric = new QualityMetric("MS:4000071", null, "number of chromatograms", t.getChromatograms().stream().count(), null);
            System.out.println("RT range...");
            var rtRange = t.getScans().stream().map(
                    scan -> Range.singleton(scan.getRetentionTime())
            ).reduce(
                    (lrt, rrt) -> lrt.span(rrt)
            ).get();
            var rtRangeMetric = new QualityMetric("MS:4000070", null, "retention time acquisition range", Arrays.asList(rtRange.lowerEndpoint(), rtRange.upperEndpoint()), new Unit(new CvParameter("UO:0000010", null, "second", null), null));
            return new SimpleEntry<InputFile, List<QualityMetric>>(
                    new InputFile(mzMLFormatParameter, Collections.emptyList(), t.getOriginalFile().get().toURI(), t.getName()),
                    Arrays.asList(
                            numberOfChromatogramsMetric,
                            ms1MzRangeMetric,
                            rtRangeMetric,
                            totalIonChromatogram
                    )
            );
        }).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

//        Plot.show(LinePlot.create("TIC plot", ticsTable, "RT [s]", "Intensity [a.u.]", "File"), new File("proteomics-qc.html"));
        try {
            var analysisSoftware = new AnalysisSoftware("MS:1000799", null, "custom unreleased software tool", "jmzqc", new URI("https://github.com/MS-Quality-hub/jmzqc"), "1.0.0-RC1");
            List<BaseQuality> bqs = mzMLFileStats.entrySet().stream().map((t) -> {
                Metadata metadata = new Metadata(Arrays.asList(analysisSoftware), Collections.emptyList(), Arrays.asList(t.getKey()), null);
                return new BaseQuality(metadata, t.getValue());
            }
            ).collect(Collectors.toList());
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
                    "MzQC for basic QC information on MetaboLights dataset MTBLS1375",
                    bqs,
                    Collections.emptyList(),
                    "1.0.0");
            Set<ValidationMessage> messages = Converter.validate(mzQC);
            System.out.println("Validation messages: " + messages);

            ObjectWriter writer = prepareJsonWriter();
            writer.writeValue(new File(outputDir, baseName + ".mzQC"), new Coordinate(mzQC));
            System.out.println(writer.writeValueAsString(new Coordinate(mzQC)).substring(0, 499));

            var ticQualityMetrics = mzQC.getRunQualityMetricsByAccession(0, "MS:1000235");
            float[] values = (float[]) ticQualityMetrics.get(0).value();
            float[] rts = (float[]) ticQualityMetrics.get(0).unit().cvParameterValue().value();
            XYChart chart = QuickChart.getChart("TIC Chart", "RT [s]", "Intensity [a.u.]",
                    mzQC.runQualities().get(0).metadata().inputFiles().get(0).name(),
                    //this plotting lib wants data as doubles, so convert our floats to double
                    IntStream.range(0, rts.length).mapToDouble(i -> rts[i]).toArray(),
                    IntStream.range(0, values.length).mapToDouble(i -> values[i]).toArray());
            //generate image
//            BitmapEncoder.getBufferedImage(chart);
            BitmapEncoder.saveBitmap(chart, new File(outputDir, "tic.png").getPath(), BitmapEncoder.BitmapFormat.PNG);

        } catch (URISyntaxException ex) {
            Logger.getLogger(JmzqcProteomics.class.getName()).log(Level.SEVERE, null, ex);
        } catch (JsonProcessingException ex) {
            Logger.getLogger(JmzqcProteomics.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(JmzqcProteomics.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    public static ObjectWriter prepareJsonWriter() {
        JsonFactoryBuilder jfb = new JsonFactoryBuilder().
                enable(JsonReadFeature.ALLOW_TRAILING_COMMA);
        ObjectMapper mapper = new ObjectMapper(jfb.build());
        mapper.findAndRegisterModules();
        mapper.configure(SerializationFeature.WRITE_DATES_AS_TIMESTAMPS, false);
        mapper.configure(SerializationFeature.INDENT_OUTPUT, true);
        mapper.configure(JsonParser.Feature.ALLOW_COMMENTS, true);
        mapper.setSerializationInclusion(Include.NON_EMPTY);

        SimpleModule module = new SimpleModule();
        module.addDeserializer(OffsetDateTime.class, new JsonDeserializer<OffsetDateTime>() {
            @Override
            public OffsetDateTime deserialize(JsonParser jsonParser, DeserializationContext deserializationContext) throws IOException, JsonProcessingException {
                String value = jsonParser.getText();
                return Converter.parseDateTimeString(value);
            }
        });
        mapper.registerModule(module);
        return mapper.writerFor(Coordinate.class);
    }

}
