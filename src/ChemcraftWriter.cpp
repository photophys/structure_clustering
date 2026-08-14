#include "ChemcraftWriter.hpp"
#include <stdexcept>

ChemcraftWriter::ChemcraftWriter(const std::string& filepath)
    : _fp(filepath)
{
    if (!_fp.is_open()) {
        throw std::runtime_error(
            "Failed to open Chemcraft output file: " + filepath
        );
    }

    // Written immediately after the file is opened.
    _fp << "Computational data saved by Chemcraft build 780\n";

    if (!_fp) {
        throw std::runtime_error(
            "Failed to write the Chemcraft file header: " + filepath
        );
    }
}

ChemcraftWriter::~ChemcraftWriter()
{
    close();
}

void ChemcraftWriter::startJob(const std::string& type, const std::string& comment)
{
    if (!_fp.is_open()) {
        throw std::runtime_error("Chemcraft output file is not open.");
    }

    if (_clusterIdx!=0)
        _fp << "[/Job]"<<std::endl;

    _clusterIdx++;
    _tempStructIdx=0;

    _fp << "[Job " << _clusterIdx << "]" << std::endl << type <<" "<<_clusterIdx;

    if (!comment.empty()) {
        _fp << " (" << comment << ')';
    }

    _fp << std::endl << "[Text]"<< std::endl
<<"[Job abstract]"<< std::endl
<<"2"<< std::endl
<<"[/Text]"<<std::endl;
}

void ChemcraftWriter::startCluster(const std::string& comment)
{
    startJob("Cluster", comment);
}

void ChemcraftWriter::startSingle(const std::string& comment)
{
    startJob("Single", comment);
}

void ChemcraftWriter::writeStructure(
    const Structure& structure,
    const std::optional<double> energy_eV)
{
    _tempStructIdx++;
    
    if (!_fp.is_open()) {
        throw std::runtime_error("Chemcraft output file is not open.");
    }



    _fp << "[Geometry]" << std::endl << "[Geometry N " << _tempStructIdx << "]"<<std::endl
        << "iso " << structure.getId();


    if (energy_eV.has_value()) {
        _fp << " (" << *energy_eV << " eV)";
    }

    _fp<<std::endl;

    _fp<<structure.numAtoms()<<std::endl;

    for (int j = 0; j < structure.numAtoms(); j++) {
        auto &atom = structure.getAtom(j);
        _fp << atom.atomicNumber() << " " << atom.position().x() << " "
             << atom.position().y() << " " << atom.position().z() << std::endl;
    }

    _fp << "[Energy: ";
    if (energy_eV.has_value()) 
        _fp << *energy_eV;
    else
        _fp << "0";
    _fp << "]"<<std::endl;

    _fp << "[/geometry]"<<std::endl<<"[Geometry]"<<std::endl;

    if (!_fp) {
        throw std::runtime_error(
            "An error occurred while writing the structure."
        );
    }
}

void ChemcraftWriter::close()
{
    if (_fp.is_open()) {
        _fp.flush();
        _fp.close();
    }
}