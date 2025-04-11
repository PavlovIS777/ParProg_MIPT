#include <iostream>
#include <cmath>
#include <mpi.h>
#include <chrono>
#include <vector>
#include <algorithm>
#include <climits>
#include <numeric>
#include <fstream>

using namespace std;

std::vector<int> read_vector_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    int size;
    file >> size;

    std::vector<int> data(size);
    for (int i = 0; i < size; ++i) {
        if (!(file >> data[i])) {
            throw std::runtime_error("Invalid data format in file");
        }
    }

    return data;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    
    auto start = chrono::steady_clock::now();
    auto end = chrono::steady_clock::now();

    int rank, pCount, dataLen;
    std::vector<int> pData;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &pCount);

    std::vector<int> delimeters;
    
    if (rank == 0)
    {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        std::vector<int> data = read_vector_from_file(argv[1]);

        start = chrono::steady_clock::now();
        std::sort(data.begin(), data.end());
        end = chrono::steady_clock::now();
        std::cout << "Elapsed time in microseconds: " << chrono::duration_cast<chrono::microseconds>(end - start).count()<< " µs" << endl;

        data = read_vector_from_file(argv[1]);

        start = chrono::steady_clock::now();
        dataLen = data.size();
        int pDataLen = data.size() / pCount;
        MPI_Bcast(&dataLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
        std::vector<MPI_Request> initSendReqs(pCount-1);

        for (int i = 1; i < pCount; i++)
        {
            MPI_Isend(
                data.data()+i*pDataLen, 
                i+1 < pCount? pDataLen : data.size()-i*pDataLen, 
                MPI_INT, 
                i, 
                i, 
                MPI_COMM_WORLD, 
                &initSendReqs[i-1]
            );
        }

        pData.resize(pDataLen);
        std::copy(data.begin(), data.begin() + pDataLen, pData.begin());
        MPI_Waitall(initSendReqs.size(), initSendReqs.data(), MPI_STATUSES_IGNORE);
    }
    else
    {
        MPI_Bcast(&dataLen, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int pDataLen = rank == pCount-1? dataLen - rank*(dataLen/pCount) : dataLen/pCount;
        pData.resize(pDataLen);
        MPI_Request dataRequest;

        MPI_Irecv(
            pData.data(), 
            pDataLen,
            MPI_INT, 
            0, 
            rank, 
            MPI_COMM_WORLD, 
            &dataRequest
        );

        MPI_Wait(&dataRequest, MPI_STATUSES_IGNORE);
    }

    std::sort(pData.begin(), pData.end());

    delimeters.resize(pCount-1);
    for (int i = 0; i < pCount; i++)
    {
        delimeters[i] = pData[pData.size()/2 + 1];
        MPI_Bcast(&delimeters[i], 1, MPI_INT, i, MPI_COMM_WORLD);
    }

    std::vector<int> sendCounts(pCount, 0);
    std::vector<int> sendOffsets(pCount, 0);

    size_t current_pos = 0;
    for (int proc = 0; proc < pCount; proc++) {
        int lower = (proc == 0) ? INT_MIN : delimeters[proc-1];
        int upper = (proc == pCount-1) ? INT_MAX : delimeters[proc];
        
        sendOffsets[proc] = current_pos;
        
        while (current_pos < pData.size() && 
            (proc == pCount-1 || pData[current_pos] <= upper)) {
            sendCounts[proc]++;
            current_pos++;
        }
    }

    std::vector<int> recvCounts(pCount);
    MPI_Alltoall(sendCounts.data(), 1, MPI_INT, recvCounts.data(), 1, MPI_INT, MPI_COMM_WORLD);
    
    int sampleLen = std::accumulate(recvCounts.begin(), recvCounts.end(), 0);

    std::vector<int> recvDispls(pCount);
    for (int i = 0; i < pCount; i++) {
        recvDispls[i] = (i == 0) ? 0 : recvDispls[i-1] + recvCounts[i-1];
    }

    std::vector<int> sample(sampleLen);
    MPI_Alltoallv(
        pData.data(), sendCounts.data(), sendOffsets.data(), MPI_INT,
        sample.data(), recvCounts.data(), recvDispls.data(), MPI_INT,
        MPI_COMM_WORLD
    );
    std::sort(sample.begin(), sample.end());

    if (rank == 0)
    {
        std::vector<int> data(dataLen);
        std::vector<int> lenghts(pCount);
        lenghts[0] = sampleLen;
        for (int i = 1; i < pCount; i++)
        {
            MPI_Recv(&lenghts[i], 1, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
        }

        std::copy(sample.begin(), sample.end(), data.begin());

        int offset = lenghts[0];
        for (int i = 1; i < pCount; i++)
        {
            MPI_Recv(data.data()+offset, lenghts[i], MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
            offset += lenghts[i];
        }

        end = chrono::steady_clock::now();
        std::cout << "Elapsed time in microseconds: " << chrono::duration_cast<chrono::microseconds>(end - start).count()<< " µs" << endl;

        std::vector<int> correctData(data.size());
        std::copy(data.begin(), data.end(), correctData.begin());

        bool valid = true;
        for(int i = 0; i < data.size(); i++)
        {
            if (data[i] != correctData[i])
            {
                std::cout << "Массивы отсортирован некорректно!";
                valid = false;
                break;
            }
        }

        if (valid)
        {
            std::cout << "Массивы отсортирован корретно!" << std::endl;
        }
        
    }
    else
    {
        MPI_Send(&sampleLen, 1, MPI_INT, 0, rank, MPI_COMM_WORLD);
        MPI_Send(sample.data(), sampleLen, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}