/* 
 * Copyright 2021 Nils Hoffmann.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.lifstools.jmzqc.usecase;

import com.fasterxml.jackson.annotation.JsonInclude;
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
import com.networknt.schema.ValidationMessage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.StringWriter;
import java.nio.file.Files;
import java.time.Instant;
import java.time.OffsetDateTime;
import java.util.Optional;
import java.util.Properties;
import java.util.Set;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.PosixParser;
import org.lifstools.jmzqc.Converter;
import org.lifstools.jmzqc.Coordinate;
import org.lifstools.jmzqc.MzQC;

/**
 * Create a new command line parser for parsing of lipid names.
 *
 * @author Dominik Kopczynski
 * @author Nils Hoffmann
 *
 */
public class CmdLineParser {

    private static String getAppInfo() throws IOException {
        Properties p = new Properties();
        p.load(CmdLineParser.class.getResourceAsStream(
                "/application.properties"));
        StringBuilder sb = new StringBuilder();
        String buildDate = p.getProperty("app.build.date", "no build date");
        if (!"no build date".equals(buildDate)) {
            Instant instant = Instant.ofEpochMilli(Long.parseLong(buildDate));
            buildDate = instant.toString();
        }
        /*
         *Property keys are in src/main/resources/application.properties
         */
        sb.append("Running ").
                append(p.getProperty("app.name", "undefined app")).
                append("\n\r").
                append(" version: '").
                append(p.getProperty("app.version", "unknown version")).
                append("'").
                append("\n\r").
                append(" build-date: '").
                append(buildDate).
                append("'").
                append("\n\r").
                append(" scm-location: '").
                append(p.getProperty("scm.location", "no scm location")).
                append("'").
                append("\n\r").
                append(" commit: '").
                append(p.getProperty("scm.commit.id", "no commit id")).
                append("'").
                append("\n\r").
                append(" branch: '").
                append(p.getProperty("scm.branch", "no branch")).
                append("'").
                append("\n\r");
        return sb.toString();
    }

    /**
     * <p>
     * Runs the command line validation for jmzqc.</p>
     *
     * Run with the {@code -h} or {@code --help} option to see more options.
     *
     * @param args an array of {@link java.lang.String} lipid names.
     * @throws java.lang.Exception if any unexpected errors occur.
     */
    @SuppressWarnings("static-access")
    public static void main(String[] args) throws Exception {
        CommandLineParser parser = new PosixParser();
        Options options = new Options();
        String helpOpt = addHelpOption(options);
        String versionOpt = addVersionOption(options);
        String inputFileOpt = addFileInputOption(options);
        String outputToFileOpt = addOutputToFileOption(options);

        CommandLine line = parser.parse(options, args);
        if (line.getOptions().length == 0 || line.hasOption(helpOpt)) {
            HelpFormatter formatter = new HelpFormatter();
            formatter.printHelp("jmzqc-usecase", options);
        } else if (line.hasOption(versionOpt)) {
            System.out.println(getAppInfo());
        } else {
            boolean toFile = false;
            String outputFile = "jmzqc.mzQC";
            if (line.hasOption(outputToFileOpt)) {
                toFile = true;
                outputFile = line.getOptionValue(outputToFileOpt);
            }
            Optional<File> inputFile = Optional.empty();
            if (line.hasOption(inputFileOpt)) {
                inputFile = Optional.of(new File(line.getOptionValue(inputFileOpt)));
            }
            Optional<MzQC> mzQC = Optional.empty();
            if (inputFile.isPresent()) {
                mzQC = new ProteomicsDDAMs1QC(inputFile.get()).process();
            }
            if (mzQC.isPresent()) {
                Set<ValidationMessage> messages = Converter.validate(mzQC.get());
                if (!messages.isEmpty()) {
                    System.out.println("Validation failed with " + messages.size() + " messages!");
                    System.out.println(messages);
                    System.exit(1);
                } else {
                    System.out.println("Validation successful!");
                }
                if (toFile) {
                    System.out.println("Saving output to '" + outputFile + "'.");
                    boolean successful = writeToFile(new File(outputFile), mzQC.get());
                    System.exit(1);
                } else {
                    System.out.println("Echoing output to stderr.");
                    boolean successful = writeToStdOut(mzQC.get());
                    System.exit(1);
                }
            } else {
                System.out.println("MzQC creation failed. Please check further output!");
                System.exit(1);
            }
        }
    }

    private static boolean writeToStdOut(MzQC mzQC) {

        try ( StringWriter sw = new StringWriter()) {
            try ( BufferedWriter bw = new BufferedWriter(sw)) {
                writeToWriter(bw, mzQC);
            }
            sw.flush();
            sw.close();
            System.err.println(sw.toString());
            return true;
        } catch (IOException ex) {
            System.err.println("Caught exception while trying to write validation results string!");
            ex.printStackTrace(System.err);
            return false;
        }
    }

    private static boolean writeToFile(File f, MzQC mzQC) {

        try ( BufferedWriter bw = Files.newBufferedWriter(f.toPath())) {
            writeToWriter(bw, mzQC);
            return true;
        } catch (IOException ex) {
            System.err.println("Caught exception while trying to write validation results to file " + f);
            ex.printStackTrace(System.err);
            return false;
        }
    }

    private static void writeToWriter(BufferedWriter bw, MzQC mzQC) {
        try {
            ObjectWriter writer = prepareJsonWriter();
            bw.write(writer.writeValueAsString(new Coordinate(mzQC)));
            bw.newLine();
        } catch (IOException ex) {
            System.err.println("Caught exception while trying to write validation results to buffered writer.");
            ex.printStackTrace(System.err);
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
        mapper.setSerializationInclusion(JsonInclude.Include.NON_EMPTY);

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

    protected static String addFileInputOption(Options options) {
        String fileOpt = "file";
        options.addOption("f", fileOpt, true, "Input file name to process.");
        return fileOpt;
    }

    protected static String addVersionOption(Options options) {
        String versionOpt = "version";
        options.addOption("v", versionOpt, false, "Print version information.");
        return versionOpt;
    }

    protected static String addHelpOption(Options options) {
        String helpOpt = "help";
        options.addOption("h", helpOpt, false, "Print help message.");
        return helpOpt;
    }

    protected static String addOutputToFileOption(Options options) {
        String outputToFileOpt = "outputFile";
        options.addOption("o", outputToFileOpt, true, "Write output to provided file in instead of to std out.");
        return outputToFileOpt;
    }

}
