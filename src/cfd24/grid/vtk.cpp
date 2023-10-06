#include "vtk.hpp"
#include <sstream>
#include <iomanip>
#include "ghc/filesystem.hpp"

using namespace cfd;

namespace fs = ghc::filesystem;

namespace{

std::vector<size_t> find_strings_in_stream(std::string fname,
                                           const std::vector<std::string>& str){
	std::ifstream fs(fname);
	std::string out(std::istreambuf_iterator<char>{fs}, std::istreambuf_iterator<char>{});
	std::vector<size_t> ret;
	for (size_t i = 0; i < str.size(); ++i)
		ret.push_back(out.find(str[i]));
	return ret;
}


std::ostream& vtk_string(std::ostream& os, const Point& p){
	os << p[0] << " " << p[1] << " " << p[2];
	return os;
}

std::ostream& vtk_string(std::ostream& s, double v){
	s << v;
	return s;
} 

std::ostream& vtk_string(std::ostream& s, const std::vector<double>& v){
	for (double d: v){
		s << d << std::endl;
	}
	return s;
}

std::ostream& vtk_string(std::ostream& s, const std::vector<Vector>& v){
	for (Vector d: v){
		s << d[0] << " " << d[1] << " " << d[2] << std::endl;
	}
	return s;
}

void create_directory(std::string path, bool purge){
	if (fs::is_directory(path)){
		if (purge){
			fs::remove_all(path);
		} else {
			return;
		}
	}
	fs::create_directory(path);
}

}

void VtkUtils::append_header(std::string caption, std::ostream& fs){
	fs << "# vtk DataFile Version 3.0" << std::endl;
	fs << caption << std::endl;
	fs << "ASCII" << std::endl;

}

void VtkUtils::append_points(const std::vector<Point>& points, std::ostream& fs){
	fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fs << "POINTS " << points.size() << " double" << std::endl;
	for (const auto& point: points){
		fs << point[0] << " " << point[1] << " " << point[2];
		fs << std::endl;
	}

}

void VtkUtils::add_cell_data(const std::vector<double>& data,
                                 std::string data_cap,
                                 std::string fname){
	// find CELL_DATA, POINT_DATA which already present
	std::vector<size_t> fnd = find_strings_in_stream(fname, {"CELL_DATA", "POINT_DATA"});

	std::ostringstream fs;
	// cell data do not exist write caption
	if (fnd[0] == std::string::npos) append_cell_data_header(data.size(), fs);

	// write data
	add_cell_data(data, data_cap, fs);

	if (fnd[1] != std::string::npos){
		// if point data exist write data before it
		std::ifstream f1(fname);
		std::string out(std::istreambuf_iterator<char>{f1}, std::istreambuf_iterator<char>{});
		f1.close();
		std::ofstream f2(fname);
		f2 << out.substr(0, fnd[1]) << fs.str() << out.substr(fnd[1], out.npos);
	} else {
		// if not, simply append
		std::ofstream f2(fname, std::ios::app);
		f2 << fs.str();
	}
}

void VtkUtils::add_cell_data(const std::vector<double>& data,
                             std::string data_cap,
                             std::ostream& fs){
	fs << "SCALARS " << data_cap << " double 1" << std::endl;
	fs << "LOOKUP_TABLE default" << std::endl;
	vtk_string(fs, data);
}

void VtkUtils::add_point_data(const std::vector<double>& data,
                              std::string data_cap,
                              std::string fname){
	// find if there are any point data already in file
	std::vector<size_t> fnd = find_strings_in_stream(fname, {"POINT_DATA"});

	// vertex data do not exist write caption
	std::fstream fs(fname, std::ios::app);
	if (fnd[0] == std::string::npos) append_point_data_header(data.size(), fs);

	// write data
	add_point_data(data, data_cap, fs);

	fs.close();
}

void VtkUtils::add_point_data(const std::vector<double>& data,
                              std::string data_cap,
                              std::ostream& fs){
	fs << "SCALARS " << data_cap << " double 1" << std::endl;
	fs << "LOOKUP_TABLE default" << std::endl;
	vtk_string(fs, data);
}
void VtkUtils::add_point_vector(const std::vector<Vector>& data,
                                std::string data_cap,
                                std::string fname){
	// find if there are any point data already in file
	std::vector<size_t> fnd = find_strings_in_stream(fname, {"POINT_DATA"});

	// vertex data do not exist write caption
	std::fstream fs(fname, std::ios::app);
	if (fnd[0] == std::string::npos) append_point_data_header(data.size(), fs);

	// write data
	add_point_vector(data, data_cap, fs);

	fs.close();
}

void VtkUtils::add_point_vector(const std::vector<Vector>& data,
                                std::string data_cap,
                                std::ostream& fs){
	fs << "SCALARS " << data_cap << " double 3" << std::endl;
	fs << "LOOKUP_TABLE default" << std::endl;
	vtk_string(fs, data);
}


void VtkUtils::append_cell_data_header(size_t data_size, std::ostream& os){
	os << "CELL_DATA " << data_size << std::endl;
}

void VtkUtils::append_point_data_header(size_t data_size, std::ostream& os){
	os << "POINT_DATA " << data_size << std::endl;
}

VtkUtils::TimeDependentWriter::TimeDependentWriter(std::string stem): _stem(stem), _series_fn(stem+".vtk.series"){
	create_directory(stem, true);

	std::ofstream ofs(_series_fn);
	if (!ofs) throw std::runtime_error("Failed to open " + _series_fn + " for writing");

	ofs << "{" << std::endl;
	ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
	ofs << "  \"files\" : [" << std::endl;
	ofs << "  ]" << std::endl;
	ofs << "}" << std::endl;
}

std::string VtkUtils::TimeDependentWriter::add(double tm){
	std::ostringstream fn;
	fn << std::setfill('0') << std::setw(8) << std::fixed << std::setprecision(4) << tm << ".vtk";
	std::string ret = _stem + '/' + fn.str();

	std::ostringstream oss;
	if (!_fileslist.empty()){
		oss << "," << std::endl;
	}
	oss << "    {\"name\": \"" << ret << "\", \"time\": " << tm << "}";
	_fileslist += oss.str();

	save_series();

	return ret;
}

void VtkUtils::TimeDependentWriter::save_series() const{
	std::fstream ofs(_series_fn);
	if (!ofs) throw std::runtime_error("Failed to open " + _series_fn + " for writing");

	ofs << "{" << std::endl;
	ofs << "  \"file-series-version\" : \"1.0\"," << std::endl;
	ofs << "  \"files\" : [" << std::endl;
	ofs << _fileslist << std::endl;
	ofs << "  ]" << std::endl;
	ofs << "}" << std::endl;
}
