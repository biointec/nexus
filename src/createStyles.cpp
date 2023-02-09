/******************************************************************************
 *  Nexus: Pan-genome compacted de Bruijn graphs with support for approximate *
 *      pattern matching using search schemes                                 *
 *                                                                            *
 *  Copyright (C) 2022 - Lore Depuydt <lore.depuydt@ugent.be>,                *
 *                       Luca Renders <luca.renders@ugent.be> and             *
 *                       Jan Fostier <jan.fostier@ugent.be>                   *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

using namespace std;

void showUsage() {
    cout << "This program creates a Cytoscape styles file that can be used to "
            "visualize subgraphs. The result wil be written to the file "
            "PanGenomeSubgraph.xml in this directory.\n\n";

    cout << "Usage: ./createStyles <numberOfStrains>\n\n";
}

int main(int argc, char* argv[]) {
    int requiredArguments = 1; // number of strains

    if (argc == 2) {
        string firstArg(argv[1]);
        if (firstArg.find("help") != std::string::npos) {
            showUsage();
            return EXIT_SUCCESS;
        }
    }

    if (argc != requiredArguments + 1) {
        std::cerr << "Error: one argument is required: the number of strains "
                     "in the pan-genome.\n"
                  << std::endl;
        showUsage();
        return EXIT_FAILURE;
    }
    std::string parameter = argv[1];

    uint32_t nr_of_strains = 0;

    try {
        nr_of_strains = std::stoi(parameter);
    } catch (const std::exception& e) {
        std::cerr
            << "Error: the parameter you provided is not a valid integer.\n"
            << std::endl;
        showUsage();
        return EXIT_FAILURE;
    }

    if (nr_of_strains == 0) {
        std::cerr << "Error: the number of strains can not be 0.\n"
                  << std::endl;
        showUsage();
        return EXIT_FAILURE;
    }

    nr_of_strains--;

    std::ofstream file;
    file.open("PanGenomeSubgraph.xml");
    file
        << "<?xml version=\"1.0\" encoding=\"UTF-8\" "
           "standalone=\"yes\"?>\n<vizmap id=\"VizMap-2022_04_22-16_47\" "
           "documentVersion=\"3.1\">\n\t<visualStyle "
           "name=\"PanGenomeSubgraph\">\n\t\t<network>\n\t\t\t<visualProperty "
           "default=\"550.0\" name=\"NETWORK_WIDTH\"/>\n\t\t\t<visualProperty "
           "default=\"\" name=\"NETWORK_TITLE\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" "
           "name=\"NETWORK_CENTER_Y_LOCATION\"/>\n\t\t\t<visualProperty "
           "default=\"false\" "
           "name=\"NETWORK_FORCE_HIGH_DETAIL\"/>\n\t\t\t<visualProperty "
           "default=\"true\" "
           "name=\"NETWORK_NODE_SELECTION\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" name=\"NETWORK_DEPTH\"/>\n\t\t\t<visualProperty "
           "default=\"400.0\" name=\"NETWORK_HEIGHT\"/>\n\t\t\t<visualProperty "
           "default=\"false\" "
           "name=\"NETWORK_ANNOTATION_SELECTION\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" "
           "name=\"NETWORK_CENTER_Z_LOCATION\"/>\n\t\t\t<visualProperty "
           "default=\"#FFFFFF\" "
           "name=\"NETWORK_BACKGROUND_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"1.0\" "
           "name=\"NETWORK_SCALE_FACTOR\"/>\n\t\t\t<visualProperty "
           "default=\"true\" "
           "name=\"NETWORK_EDGE_SELECTION\"/>\n\t\t\t<visualProperty "
           "default=\"false\" "
           "name=\"NETWORK_NODE_LABEL_SELECTION\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" "
           "name=\"NETWORK_CENTER_X_LOCATION\"/>\n\t\t\t<visualProperty "
           "default=\"550.0\" "
           "name=\"NETWORK_SIZE\"/>\n\t\t</"
           "network>\n\t\t<node>\n\t\t\t<dependency value=\"true\" "
           "name=\"nodeCustomGraphicsSizeSync\"/>\n\t\t\t<dependency "
           "value=\"false\" name=\"nodeSizeLocked\"/>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_2, "
           "name=Node Custom Paint 2)\" "
           "name=\"NODE_CUSTOMPAINT_2\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" name=\"NODE_DEPTH\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" "
           "name=\"NODE_BORDER_WIDTH\"/>\n\t\t\t<visualProperty "
           "default=\"#000000\" "
           "name=\"NODE_LABEL_COLOR\"/>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_7, "
           "name=Node Custom Paint 7)\" "
           "name=\"NODE_CUSTOMPAINT_7\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_1\"/>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_9, "
           "name=Node Custom Paint 9)\" "
           "name=\"NODE_CUSTOMPAINT_9\"/>\n\t\t\t<visualProperty "
           "default=\"ROUND_RECTANGLE\" "
           "name=\"NODE_SHAPE\"/>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_8, "
           "name=Node Custom Paint 8)\" "
           "name=\"NODE_CUSTOMPAINT_8\"/>\n\t\t\t<visualProperty "
           "default=\"false\" name=\"NODE_SELECTED\"/>\n\t\t\t<visualProperty "
           "default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_3\"/>\n\t\t\t<visualProperty "
           "default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_5\"/>\n\t\t\t<visualProperty "
           "default=\"#1E90FF\" name=\"NODE_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_3\"/>\n\t\t\t<visualProperty "
           "default=\"200.0\" "
           "name=\"NODE_LABEL_WIDTH\"/>\n\t\t\t<visualProperty "
           "default=\"35.0\" name=\"NODE_SIZE\"/>\n\t\t\t<visualProperty "
           "default=\"\" name=\"NODE_LABEL\">\n\t\t\t\t<passthroughMapping "
           "attributeName=\"OmegaShort\" "
           "attributeType=\"string\"/>\n\t\t\t</"
           "visualProperty>\n\t\t\t<visualProperty default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_8\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_6\"/>\n\t\t\t<visualProperty "
           "default=\"SOLID\" "
           "name=\"NODE_BORDER_STROKE\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_4\"/>\n\t\t\t<visualProperty "
           "default=\"#FFFF00\" "
           "name=\"NODE_SELECTED_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_6\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" name=\"NODE_Y_LOCATION\"/>\n\t\t\t<visualProperty "
           "default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_9\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_6\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_7\"/>\n\t\t\t<visualProperty "
           "default=\"true\" "
           "name=\"NODE_NESTED_NETWORK_IMAGE_VISIBLE\"/"
           ">\n\t\t\t<visualProperty default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_1\"/>\n\t\t\t<visualProperty "
           "default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_4\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_LABEL_POSITION\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_8\"/>\n\t\t\t<visualProperty "
           "default=\"255\" "
           "name=\"NODE_LABEL_TRANSPARENCY\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_7\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_9\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" name=\"NODE_Z_LOCATION\"/>\n\t\t\t<visualProperty "
           "default=\"12\" "
           "name=\"NODE_LABEL_FONT_SIZE\"/>\n\t\t\t<visualProperty "
           "default=\"\" name=\"NODE_TOOLTIP\"/>\n\t\t\t<visualProperty "
           "default=\"#000000\" "
           "name=\"NODE_BORDER_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"255\" "
           "name=\"NODE_BORDER_TRANSPARENCY\"/>\n\t\t\t<visualProperty "
           "default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_2\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_5\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_9\"/>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_5, "
           "name=Node Custom Paint 5)\" "
           "name=\"NODE_CUSTOMPAINT_5\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_2\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" name=\"NODE_X_LOCATION\"/>\n\t\t\t<visualProperty "
           "default=\"255\" "
           "name=\"NODE_TRANSPARENCY\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" "
           "name=\"NODE_LABEL_ROTATION\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_5\"/>\n\t\t\t<visualProperty "
           "default=\"50.0\" "
           "name=\"NODE_CUSTOMGRAPHICS_SIZE_7\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_8\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_3\"/>\n\t\t\t<visualProperty "
           "default=\"35.0\" name=\"NODE_HEIGHT\"/>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_1, "
           "name=Node Custom Paint 1)\" "
           "name=\"NODE_CUSTOMPAINT_1\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_4\"/>\n\t\t\t<visualProperty "
           "default=\"ROUND_RECTANGLE\" "
           "name=\"COMPOUND_NODE_SHAPE\"/>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_6, "
           "name=Node Custom Paint 6)\" "
           "name=\"NODE_CUSTOMPAINT_6\"/>\n\t\t\t<visualProperty "
           "default=\"10.0\" "
           "name=\"COMPOUND_NODE_PADDING\"/>\n\t\t\t<visualProperty "
           "default=\"#89D0F5\" "
           "name=\"NODE_FILL_COLOR\">\n\t\t\t\t<discreteMapping "
           "attributeName=\"PartOfPath\" "
           "attributeType=\"integer\">\n\t\t\t\t\t<discreteMappingEntry "
           "attributeValue=\"0\" "
           "value=\"#D8D8D8\"/>\n\t\t\t\t\t<discreteMappingEntry "
           "attributeValue=\"1\" "
           "value=\"#A1A1A1\"/>\n\t\t\t\t</discreteMapping>\n\t\t\t</"
           "visualProperty>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_4, "
           "name=Node Custom Paint 4)\" "
           "name=\"NODE_CUSTOMPAINT_4\"/>\n\t\t\t<visualProperty "
           "default=\"SansSerif.plain,plain,12\" "
           "name=\"NODE_LABEL_FONT_FACE\"/>\n\t\t\t<visualProperty "
           "default=\"DefaultVisualizableVisualProperty(id=NODE_CUSTOMPAINT_3, "
           "name=Node Custom Paint 3)\" "
           "name=\"NODE_CUSTOMPAINT_3\"/>\n\t\t\t<visualProperty "
           "default=\"org.cytoscape.cg.model.NullCustomGraphics,0,[ Remove "
           "Graphics ],\" "
           "name=\"NODE_CUSTOMGRAPHICS_1\"/>\n\t\t\t<visualProperty "
           "default=\"C,C,c,0.00,0.00\" "
           "name=\"NODE_CUSTOMGRAPHICS_POSITION_2\"/>\n\t\t\t<visualProperty "
           "default=\"true\" name=\"NODE_VISIBLE\"/>\n\t\t\t<visualProperty "
           "default=\"250.0\" "
           "name=\"NODE_WIDTH\"/>\n\t\t</node>\n\t\t<edge>\n\t\t\t<dependency "
           "value=\"true\" "
           "name=\"arrowColorMatchesEdge\"/>\n\t\t\t<visualProperty "
           "default=\"6.0\" "
           "name=\"EDGE_TARGET_ARROW_SIZE\"/>\n\t\t\t<visualProperty "
           "default=\"\" name=\"EDGE_TOOLTIP\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" "
           "name=\"EDGE_LABEL_ROTATION\"/>\n\t\t\t<visualProperty "
           "default=\"255\" "
           "name=\"EDGE_LABEL_TRANSPARENCY\"/>\n\t\t\t<visualProperty "
           "default=\"#FF0000\" "
           "name=\"EDGE_SELECTED_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"#848484\" "
           "name=\"EDGE_STROKE_UNSELECTED_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"#323232\" name=\"EDGE_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"NONE\" "
           "name=\"EDGE_SOURCE_ARROW_SHAPE\"/>\n\t\t\t<visualProperty "
           "default=\"true\" name=\"EDGE_CURVED\"/>\n\t\t\t<visualProperty "
           "default=\"SOLID\" name=\"EDGE_LINE_TYPE\"/>\n\t\t\t<visualProperty "
           "default=\"2.0\" name=\"EDGE_WIDTH\"/>\n\t\t\t<visualProperty "
           "default=\"0.0\" name=\"EDGE_Z_ORDER\"/>\n\t\t\t<visualProperty "
           "default=\"\" name=\"EDGE_LABEL\"/>\n\t\t\t<visualProperty "
           "default=\"10\" "
           "name=\"EDGE_LABEL_FONT_SIZE\"/>\n\t\t\t<visualProperty "
           "default=\"255\" "
           "name=\"EDGE_TRANSPARENCY\"/>\n\t\t\t<visualProperty "
           "default=\"#000000\" "
           "name=\"EDGE_TARGET_ARROW_UNSELECTED_PAINT\"/"
           ">\n\t\t\t<visualProperty default=\"Dialog.plain,plain,10\" "
           "name=\"EDGE_LABEL_FONT_FACE\"/>\n\t\t\t<visualProperty "
           "default=\"\" name=\"EDGE_BEND\"/>\n\t\t\t<visualProperty "
           "default=\"0.5\" "
           "name=\"EDGE_STACKING_DENSITY\"/>\n\t\t\t<visualProperty "
           "default=\"#FFFF00\" "
           "name=\"EDGE_TARGET_ARROW_SELECTED_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"#000000\" "
           "name=\"EDGE_SOURCE_ARROW_UNSELECTED_PAINT\"/"
           ">\n\t\t\t<visualProperty default=\"#FF0000\" "
           "name=\"EDGE_STROKE_SELECTED_PAINT\"/>\n\t\t\t<visualProperty "
           "default=\"#000000\" "
           "name=\"EDGE_LABEL_COLOR\"/>\n\t\t\t<visualProperty "
           "default=\"ARROW\" "
           "name=\"EDGE_TARGET_ARROW_SHAPE\">\n\t\t\t\t<continuousMapping "
           "attributeName=\"Color\" "
           "attributeType=\"float\"/>\n\t\t\t</"
           "visualProperty>\n\t\t\t<visualProperty default=\"200.0\" "
           "name=\"EDGE_LABEL_WIDTH\"/>\n\t\t\t<visualProperty "
           "default=\"AUTO_BEND\" "
           "name=\"EDGE_STACKING\"/>\n\t\t\t<visualProperty default=\"6.0\" "
           "name=\"EDGE_SOURCE_ARROW_SIZE\"/>\n\t\t\t<visualProperty "
           "default=\"true\" name=\"EDGE_VISIBLE\"/>\n\t\t\t<visualProperty "
           "default=\"false\" name=\"EDGE_SELECTED\"/>\n\t\t\t<visualProperty "
           "default=\"#404040\" "
           "name=\"EDGE_UNSELECTED_PAINT\">\n\t\t\t\t<continuousMapping "
           "attributeName=\"Color\" "
           "attributeType=\"float\">\n\t\t\t\t\t<continuousMappingPoint "
           "attrValue=\"0.0\" equalValue=\"#FF00FF\" greaterValue=\"#FF00FF\" "
           "lesserValue=\"#000000\"/>\n\t\t\t\t\t<continuousMappingPoint "
           "attrValue=\"";
    file << std::fixed << std::setprecision(1) << ((float)nr_of_strains / 2.0);
    file << "\" equalValue=\"#FFFF00\" greaterValue=\"#FFFF00\" "
            "lesserValue=\"#FFFF00\"/>\n\t\t\t\t\t<continuousMappingPoint "
            "attrValue=\"";
    file << std::fixed << std::setprecision(1) << (float)nr_of_strains;
    file << "\" equalValue=\"#00FFFF\" greaterValue=\"#000000\" "
            "lesserValue=\"#00FFFF\"/>\n\t\t\t\t</continuousMapping>\n\t\t\t</"
            "visualProperty>\n\t\t\t<visualProperty default=\"#FFFF00\" "
            "name=\"EDGE_SOURCE_ARROW_SELECTED_PAINT\"/>\n\t\t</edge>\n\t</"
            "visualStyle>\n</vizmap>\n";

    file.close();

    cout << "Program has terminated successfully" << endl;
}
