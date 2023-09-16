#ifndef CFD24_GEOM_VTK_HPP
#define CFD24_GEOM_VTK_HPP

#include "cfd24/grid/i_grid.hpp"
#include <fstream>

namespace cfd{

/**
 * @brief Collection of vtk-writer utilities 
 */
struct VtkUtils{
	/// adds vtk caption
	static void append_header(std::string caption, std::ostream& fs);

	/// adds point list
	static void append_points(const std::vector<Point>& points, std::ostream& fs);

	/**
	 * @brief adds cell data to saved vtk grid
	 * @param data data vector
	 * @param data_cap data caption as it should appear in vtk
	 * @param fname file name of saved vtk grid
	 */
	static void add_cell_data(const std::vector<double>& data, std::string data_cap, std::string fname);

	/**
	 * @brief adds cell data to a stream
	 *
	 * @param data data vector
	 * @param data_cap data caption as it should appear in vtk
	 * @param s an opened stream
	 *
	 * @note this method do not paste header of cell data section
	 * and don't care that insertion will be done into correct place
	 */
	static void add_cell_data(const std::vector<double>& data, std::string data_cap, std::ostream& s);

	/**
	 * @brief froms and saves header of cell data section
	 *
	 * @param data_size data size which equals cells number
	 * @param s an opened stream
	 */
	static void append_cell_data_header(size_t data_size, std::ostream& s);

	/**
	 * @brief adds vertex data to saved vtk file
	 * @param data data vector
	 * @param data_cap data caption as it should appear in vtk
	 * @param fname saved vtk file name
	 */
	static void add_point_data(const std::vector<double>& data,
	                           std::string data_cap,
	                           std::string fname);

	/**
	 * @brief adds vertex data to a stream
	 *
	 * @param data data vector
	 * @param data_cap data caption as it should appear in vtk
	 * @param s an opened stream
	 *
	 * @note this method do not paste header of point data section
	 *       and don't care that insertion will be done into correct place
	 */
	static void add_point_data(const std::vector<double>& data, std::string data_cap, std::ostream& s);

	/**
	 * @brief froms and saves header of point data section
	 *
	 * @param data_size data size which equals points number
	 * @param s an opened stream
	 */
	static void append_point_data_header(size_t data_size, std::ostream& s);

	class TimeDependentWriter{
	public:
		TimeDependentWriter(std::string stem);
		std::string add(double tm);
	private:
		const std::string _stem;
		const std::string _series_fn;
		bool _first_entry = true;
	};
};

}

#endif
