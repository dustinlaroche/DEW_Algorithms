//
// Created by dlaroche on 6/7/24.
//

#ifndef SIGNAL_SOURCING_H
#define SIGNAL_SOURCING_H
#include <iostream>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath> // for sin, cos, degrees
#endif //SIGNAL_SOURCING_H

#define EARTH_RADIUS_FEET (3963 * std::pow(10, 3)) // Earth radius in feet (approx)


typedef complex<double> Complex;

// Earth's equatorial radius in meters (WGS84)
const double equatorial_radius = 6378137.0;
const double EARTH_RADIUS_METERS = 6371000.0;


struct SignalInfo {
    size_t index_frequency;
    double wavelength;
    double horizontal_angle_degrees;
    double horizontal_angle_magnetic_degrees;
    double vertical_angle_degrees;
    double source_altitude;
    double distance_imag_feet;
    double distance_real_feet;
    double distance_imag_meters;
    double distance_real_meters;
    double horizontal_radians;
};



// Function to get the magnitudes of the FFT results
std::vector<double> get_magnitudes(std::vector<Complex> fft_result, double adc_gain) {
    std::vector<double> magnitudes(fft_result.size());

    for (int i = 0; i < (int)fft_result.size(); i++) {
        magnitudes[i] = sqrt(real(fft_result[i]) * real(fft_result[i]) + imag(fft_result[i]) * imag(fft_result[i])) / adc_gain;
    }

    return magnitudes;
}


// Function to find frequencies above a given threshold
std::vector<double> find_frequencies_above_threshold(std::vector<double> magnitudes, double threshold) {
    std::vector<double> frequencies;
    for (int i = 0; i < (int)magnitudes.size(); i++) {
        if (magnitudes[i] > threshold) {
            frequencies.push_back(i);
        }
    }
    return frequencies;
}


// Function to convert frequency to wavelength
double calculateWavelength(double frequency) {
    const double speedOfLight = 299792458.0; // meters per second
    return speedOfLight / frequency;
}


double DOA_Horizontal_Angle(std::complex<double> fft_result, bool fromMagneticNorth = false) {
    double angle_rad = std::arg(fft_result);

    double angle_deg = angle_rad * 180.0 / M_PI;

    if (fromMagneticNorth) {
        angle_deg = fmod(angle_deg + 360.0, 360.0);
    }

    return angle_deg;
}


double calculateAltitude(double receivingAltitude, double distanceToOrigin,  double originVerticalAngle, bool seaLevel = true) {
    if(seaLevel == false)
        return receivingAltitude + distanceToOrigin * tan(originVerticalAngle * M_PI / 180);
    else
        return distanceToOrigin * tan(originVerticalAngle * M_PI / 180) - receivingAltitude;
}


double calculateElevation(const std::vector<std::complex<double>>& fft_result, double magnitude_threshold) {
    for (const auto& sample : fft_result) {
        double magnitude = std::abs(sample);

        if (magnitude > magnitude_threshold) {
            // Assuming the elevation is related to the phase of the complex sample
            double elevation = std::atan2(std::imag(sample), std::real(sample));

            // Return the elevation in degrees
            return elevation;
        }
    }

    // If no samples are above the threshold, return a default value (e.g., 0)
    return 0.0;
}


double calculateDistance(double magnitude, double wavelength) {
    const double speedOfLight = 299792458.0; // meters per second
    double frequency = speedOfLight / wavelength; // calculate frequency
    return magnitude / (2 * M_PI * frequency) * 3.28084; // calculate distance
}


/*
// Function to calculate the signal's geo-coordinates
void calculateGeoCoordinates(const SignalInfo& signal, double lat_rad, double lon_rad, double& out_lat_rad, double& out_lon_rad) {
    // Convert angles to radians
    double horizontal_rad = deg2rad(signal.horizontal_angle_degrees);
    double vertical_rad = deg2rad(signal.vertical_angle_degrees);

    // Calculate distance along horizontal plane
    double horizontal_distance = signal.distance_feet * cos(vertical_rad);

    // Calculate offsets in latitude and longitude based on distance and angles
    double lat_offset = horizontal_distance * sin(horizontal_rad) / EARTH_RADIUS_FEET;
    double lon_offset = horizontal_distance * cos(horizontal_rad) / (EARTH_RADIUS_FEET * cos(lat_rad));

    // Add offsets to receiving coordinates
    out_lat_rad = lat_rad + lat_offset;
    out_lon_rad = lon_rad + lon_offset;
}
*/

/*
// Function to calculate the origin point coordinates
std::pair<double, double> calculate_origin(double receiver_latitude, double receiver_longitude, double horizontal_angle_deg, double distance_meters) {
    // Convert angle to radians
    double horizontal_angle_rad = deg2rad(horizontal_angle_deg);

    // Calculate initial position based on distance and angle
    double origin_latitude = receiver_latitude + distance_meters * cos(horizontal_angle_rad) / equatorial_radius;
    double origin_longitude = receiver_longitude + distance_meters * sin(horizontal_angle_rad) * cos(deg2rad(receiver_latitude)) / equatorial_radius;

    // Iterate to account for Earth's curvature (optional for short distances)
    // This loop can be omitted for most practical cases
    // for (int i = 0; i < 10; ++i) {
        double delta_latitude = origin_latitude - receiver_latitude;
        double delta_longitude = origin_longitude - receiver_longitude;
        double updated_latitude = receiver_latitude + distance_meters * cos(horizontal_angle_rad) / sqrt(equatorial_radius * equatorial_radius + delta_latitude * delta_latitude * cos(deg2rad(receiver_latitude)) * cos(deg2rad(receiver_latitude)));
        double updated_longitude = receiver_longitude + distance_meters * sin(horizontal_angle_rad) * cos(deg2rad(receiver_latitude)) / sqrt(equatorial_radius * equatorial_radius + delta_longitude * delta_longitude * cos(deg2rad(receiver_latitude)) * cos(deg2rad(receiver_latitude)));
        origin_latitude = updated_latitude;
        origin_longitude = updated_longitude;
   // }

    return std::make_pair(origin_latitude, origin_longitude);
}
*/

/*
std::pair<double, double> next_geo_coordinate(double lat, double lon, double bearing, double distance) {
    double lat1 = deg2rad(lat);
    double lon1 = deg2rad(lon);
    double distance_rad = distance / equatorial_radius;

    double lat2 = asin(sin(lat1)) * cos(distance_rad) + cos(lat1) * sin(distance_rad) * cos(bearing);
    double lon2 = lon1 + atan2(sin(bearing) * sin(distance_rad) * cos(lat1), cos(distance_rad) - sin(lat1) * sin(lat2));

    double origin_lat = rad2deg(lat2);
    double origin_lon = rad2deg(lon2);

    return std::make_pair(origin_lat, origin_lon);
}
*/

/*
std::pair<double, double> next_geo_coordinate2(double lat, double lon, double bearing, double distance) {
    double lat1 = deg2rad(lat);
    double lon1 = deg2rad(lon);
    double distance_rad = distance / equatorial_radius;

    double lat2 = asin(sin(lat1) * cos(distance_rad) + cos(lat1) * sin(distance_rad) * cos(bearing));
    double lon2 = lon1 + atan2(sin(bearing) * sin(distance_rad) * cos(lat1), cos(distance_rad) - sin(lat1) * sin(lat2));

    double origin_lat = rad2deg(lat2);
    double origin_lon = rad2deg(lon2);

    return std::make_pair(origin_lat, origin_lon);
}
*/

/*
void createKML(string file_path_name, double receiver_latitude_deg, double receiver_longitude_deg, vector<SignalInfo> signals) {
    // Open KML file for writing
    std::ofstream kml_file(file_path_name);
    if (!kml_file.is_open()) {
        std::cerr << "Error opening KML file!" << std::endl;
        return;
    }

    // Write KML header
    kml_file << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
             << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n"
             << "<Document>\n";

    // Loop through signals and write placemarks
    for (SignalInfo& signal : signals) {
        // std::pair<double, double> origin_point = calculate_origin(receiver_latitude_deg, receiver_longitude_deg, signal.horizontal_angle_degrees, signal.distance_meters);
        std::pair<double, double> origin_point = next_geo_coordinate2(receiver_latitude_deg, receiver_longitude_deg, signal.horizontal_radians, signal.distance_imag_meters);

        // Write placemark for each signal
        kml_file << "<Placemark>\n"
                 << "<name>Signal " << signal.index_frequency << "</name>\n"
                 << "<description>Wavelength: " << signal.wavelength << " meters, "
                 << "Distance: " << signal.distance_imag_meters << " feet</description>\n"
                 << "<Point>\n"
                 << "<coordinates>" << origin_point.first << ","
                 << origin_point.second << "," << signal.source_altitude << "</coordinates>\n"
                 << "</Point>\n"
                 << "</Placemark>\n";
    }

    kml_file << "</Document>\n"
             << "</kml>\n";
}
*/




double DOA_Vertical_Angle(Complex fft_result) {
    return fmod(std::atan2(std::imag(fft_result), std::real(fft_result)) + 90.0, -90.0);
    //asin(phase_diff / (2 * M_PI * d / lambda)) * 180.0 / M_PI; // Convert to degrees
}


std::pair<double, double> DOA_Distance_Feet(Complex fft_result, double frequency) {
    double phaseShift = std::arg(fft_result); // Phase shift in radians
    double distance = (phaseShift * calculateWavelength(frequency)) / (2 * M_PI);
    double real_distance = abs(distance) * 3.28084; // + 100000; // convert to feet
    double imag_distance = abs(distance) * 3.28084; // + 100000; // convert to feet

    return std::make_pair(imag_distance, real_distance);
}


std::pair<double, double> DOA_Distance_Meters(Complex fft_result, double frequency) {
    double phaseShift = std::arg(fft_result.real()); // Phase shift in radians
    double real_distance = abs((phaseShift * calculateWavelength(frequency)) / (2 * M_PI));

    phaseShift = std::arg(fft_result.imag()); // Phase shift in radians
    double imag_distance = abs((phaseShift * calculateWavelength(frequency)) / (2 * M_PI));

    return std::make_pair(imag_distance, real_distance);
}


double DOA_Source_Altitude(Complex fft_result, double frequency, double receiverMetersAboveSeaLevel = 539) {
    double altitude = receiverMetersAboveSeaLevel + DOA_Distance_Meters(fft_result, frequency).first * std::tan(DOA_Vertical_Angle(fft_result) * M_PI / 180);
    return altitude;
}


const double PI = 3.141592653589793;
const double SPEED_OF_LIGHT = 299792458.0;  // Speed of light in m/s

double calculatePhaseDifference(const std::complex<double>& sample1, const std::complex<double>& sample2) {
    return std::arg(sample2) - std::arg(sample1);
}


double calculateAngleOfArrival(double phaseDifference, double wavelength, double antennaDistance) {
    return std::asin((phaseDifference * wavelength) / (2 * PI * antennaDistance));
}


// This one will work with both polarity samples with DC removed and needs antenna spacing
double calculateBearingFromNorth(const std::vector<std::complex<double>>& samplesAntenna1,
                                 const std::vector<std::complex<double>>& samplesAntenna2,
                                 double frequency, double antennaDistance) {
    if (samplesAntenna1.size() != samplesAntenna2.size()) {
        throw std::invalid_argument("Sample arrays must have the same size.");
    }

    size_t sampleCount = samplesAntenna1.size();
    double phaseDifferenceSum = 0.0;

    for (size_t i = 0; i < sampleCount; ++i) {
        phaseDifferenceSum += calculatePhaseDifference(samplesAntenna1[i], samplesAntenna2[i]);
    }

    double averagePhaseDifference = phaseDifferenceSum / sampleCount;
    double wavelength = calculateWavelength(frequency);
    double angleOfArrival = calculateAngleOfArrival(averagePhaseDifference, wavelength, antennaDistance);

    // Convert angle of arrival to bearing from magnetic north
    double bearing = angleOfArrival * (180.0 / PI);  // Convert to degrees

    // Ensure bearing is within [0, 360) degrees
    if (bearing < 0) {
        bearing += 360.0;
    }

    return bearing;
}



// Function to convert degrees to radians
double degToRad(double degrees) {
    return degrees * M_PI / 180.0;
}

// Function to convert radians to degrees
double radToDeg(double radians) {
    return radians * 180.0 / M_PI;
}

// Helper function to generate KML placemark string
std::string createPlacemark(double lat, double lon, double altitude, const std::string& name, const std::string& description) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(6);
    oss << "<Placemark>\n";
    oss << "  <name>" << name << "</name>\n";
    oss << "  <description>" << description << " Altitude: " << altitude << "</description>\n";
    oss << "  <Point>\n";
    oss << "    <coordinates>" << lon << "," << lat << "," << altitude << "</coordinates>\n";
    oss << "  </Point>\n";
    oss << "</Placemark>\n";
    return oss.str();
}

std::pair<double, double> calculateDestination(double startLat, double startLon, double distanceMeters, double bearingDegrees) {
    double bearingRad = degToRad(bearingDegrees);
    double startLatRad = degToRad(startLat);
    double startLonRad = degToRad(startLon);

    double destLatRad = asin(sin(startLatRad) * cos(distanceMeters / EARTH_RADIUS_METERS) +
                             cos(startLatRad) * sin(distanceMeters / EARTH_RADIUS_METERS) * cos(bearingRad));

    double destLonRad = startLonRad + atan2(sin(bearingRad) * sin(distanceMeters / EARTH_RADIUS_METERS) * cos(startLatRad),
                                            cos(distanceMeters / EARTH_RADIUS_METERS) - sin(startLatRad) * sin(destLatRad));

    double destLat = radToDeg(destLatRad);
    double destLon = radToDeg(destLonRad);

    return std::make_pair(destLat, destLon);
}


void createKMLFile(const std::string& filename, const std::vector<SignalInfo>& signals, double startLat, double startLon) {
    std::ofstream kmlFile(filename);
    if (!kmlFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // KML header
    kmlFile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    kmlFile << "<kml xmlns=\"http://www.opengis.net/kml/2.2\">\n";
    kmlFile << "<Document>\n";
    kmlFile << "  <name>Signals.kml</name>\n";

    // Add placemarks for each signal
    for (const auto& signal : signals) {
        double altitude = signal.source_altitude;
        std::string name = "Signal at Frequency " + std::to_string(signal.index_frequency);
        std::string description = "Horizontal Angle: " + std::to_string(signal.horizontal_angle_degrees) + " degrees, " +
                                  "Vertical Angle: " + std::to_string(signal.vertical_angle_degrees) + " degrees";

        double distance = signal.distance_real_meters; // Use the real distance for horizontal distance
        double bearing = signal.horizontal_angle_degrees; // Use the horizontal angle for bearing

        auto [destLat, destLon] = calculateDestination(startLat, startLon, distance, bearing);

        kmlFile << createPlacemark(destLat, destLon, altitude, name, description);
    }

    // KML footer
    kmlFile << "</Document>\n";
    kmlFile << "</kml>\n";

    kmlFile.close();
    std::cout << "KML file created successfully: " << filename << std::endl;
}



// Find signals above threshold in FFT results with their direction angle and distance
std::vector<SignalInfo> find_signals_above_threshold(const std::vector<std::complex<double>>& fft_result, int threshold) {
    std::vector<SignalInfo> result;

    for (size_t frequency_bin = 0; frequency_bin < fft_result.size(); frequency_bin++) {
        double magnitude = abs(std::real(fft_result[frequency_bin]));

        if (magnitude > threshold) {
            SignalInfo signal_info;
            signal_info.index_frequency = frequency_bin;
            signal_info.wavelength = calculateWavelength(frequency_bin);
            signal_info.horizontal_angle_degrees = DOA_Horizontal_Angle(fft_result[frequency_bin], false);
            signal_info.horizontal_angle_magnetic_degrees = DOA_Horizontal_Angle(fft_result[frequency_bin], true); // normalize to 0-360 range
            signal_info.vertical_angle_degrees = DOA_Vertical_Angle(fft_result[frequency_bin]);
            signal_info.distance_imag_feet = DOA_Distance_Feet(fft_result[frequency_bin], frequency_bin).first;
            signal_info.distance_real_feet = DOA_Distance_Feet(fft_result[frequency_bin], frequency_bin).second;
            signal_info.distance_imag_meters = DOA_Distance_Meters(fft_result[frequency_bin], frequency_bin).first;
            signal_info.distance_real_meters = DOA_Distance_Meters(fft_result[frequency_bin], frequency_bin).second;
            signal_info.source_altitude = DOA_Source_Altitude(fft_result[frequency_bin], frequency_bin);
            signal_info.horizontal_radians = std::arg(fft_result[frequency_bin]);
            result.push_back(signal_info);
        }
    }

    return result;
}
