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
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.time.OffsetDateTime;
import java.util.AbstractMap.SimpleEntry;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
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

/**
 *
 * @author nilshoffmann
 */
public class JmzqcUsecase {

    public static void main(String[] args) {
        var outputDir = new File("MTBLS1375");
        outputDir.mkdirs();
        try {
            var cmd = "wget -r -l 1 -A .mzML --reject 'index.html*' -P MTBLS1375 -nc -np -nd -nH --cut-dirs=6 -e robots=off https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS1375/";
            var pb = new ProcessBuilder(cmd.split(" "));
            System.out.println("Running command: " + pb.command());
            pb.inheritIO();
            Process p = pb.start();
            p.waitFor();
        } catch (IOException | InterruptedException ex) {
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

            var precursorMzRange = t.getChromatograms().stream().filter(chrom -> chrom.getChromatogramType() == ChromatogramType.MRM_SRM).map(chrom -> {
                var mzr = Range.singleton(chrom.getIsolations().get(0).getPrecursorMz());
                return mzr;
            }).filter(range -> range != null && !range.isEmpty()).reduce((l, r) -> l.span(r)).orElse(Range.singleton(0.0));
            var precursorMzRangeMetric = new QualityMetric("MS:4000069", null, "m/z acquisition range", Arrays.asList(precursorMzRange.lowerEndpoint(), precursorMzRange.upperEndpoint()), null);
            var numberOfChromatogramsMetric = new QualityMetric("MS:4000071", null, "number of chromatograms", t.getChromatograms().stream().filter(chrom -> chrom.getChromatogramType() == ChromatogramType.MRM_SRM).count(), null);
            var rtRange = t.getChromatograms().stream().map(
                    chrom -> {
                        return chrom.getRtRange();
                    }
            ).reduce(
                    (lrt, rrt) -> lrt.span(rrt)
            ).get();
            var rtRangeMetric = new QualityMetric("MS:4000070", null, "retention time acquisition range", Arrays.asList(rtRange.lowerEndpoint(), rtRange.upperEndpoint()), null);
            return new SimpleEntry<InputFile, List<QualityMetric>>(
                    new InputFile(mzMLFormatParameter, Collections.emptyList(), t.getOriginalFile().get().toURI(), t.getName()),
                    Arrays.asList(
                            numberOfChromatogramsMetric,
                            precursorMzRangeMetric,
                            rtRangeMetric
                    )
            );
        }).collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));

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
            writer.writeValue(new File("MTBLS1375.mzQC"), new Coordinate(mzQC));
        } catch (URISyntaxException ex) {
            Logger.getLogger(JmzqcUsecase.class.getName()).log(Level.SEVERE, null, ex);
        } catch (JsonProcessingException ex) {
            Logger.getLogger(JmzqcUsecase.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(JmzqcUsecase.class.getName()).log(Level.SEVERE, null, ex);
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
