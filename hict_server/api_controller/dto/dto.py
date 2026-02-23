#  MIT License
#
#  Copyright (c) 2021-2026. Aleksandr Serdiukov, Anton Zamyatin, Aleksandr Sinitsyn, Vitalii Dravgelis and Computer Technologies Laboratory ITMO University team.
#
#  Permission is hereby granted, free of charge, to any person obtaining a copy
#  of this software and associated documentation files (the "Software"), to deal
#  in the Software without restriction, including without limitation the rights
#  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#  copies of the Software, and to permit persons to whom the Software is
#  furnished to do so, subject to the following conditions:
#
#  The above copyright notice and this permission notice shall be included in all
#  copies or substantial portions of the Software.
#
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#  SOFTWARE.

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple
import numpy as np
from hict.core.common import ContigDescriptor, ScaffoldDescriptor, ScaffoldBordersBP, ContigDirection


@dataclass
class ContigDescriptorDTO:
    contigId: int
    contigName: str
    contigDirection: int
    contigLengthBp: int
    contigLengthBins: Dict[int, int]
    # scaffoldId: Optional[int]
    contigPresenceAtResolution: Dict[int, int]

    @staticmethod
    def fromEntity(descriptor: ContigDescriptor, direction: ContigDirection) -> 'ContigDescriptorDTO':
        contig_length_at_resolution: Dict[int, int] = dict()
        presence_in_resolution: Dict[int, int] = dict()

        for res, ctg_length in descriptor.contig_length_at_resolution.items():
            if res != 0:
                int_res: int = int(res)
                contig_length_at_resolution[int_res] = int(ctg_length)
                presence_in_resolution[int_res] = descriptor.presence_in_resolution[res].value

        return ContigDescriptorDTO(
            contigId=int(descriptor.contig_id),
            contigName=str(descriptor.contig_name),
            contigDirection=direction.value,
            contigLengthBp=int(descriptor.contig_length_at_resolution[0]),
            contigLengthBins=contig_length_at_resolution,
            contigPresenceAtResolution=presence_in_resolution
        )


@dataclass
class ScaffoldBordersBPDTO:
    startBP: int
    endBP: int

    @staticmethod
    def fromEntity(borders: Optional[ScaffoldBordersBP]) -> Optional['ScaffoldBordersBPDTO']:
        return ScaffoldBordersBPDTO(
            int(borders.start_bp),
            int(borders.end_bp)
        ) if borders is not None else None


@dataclass
class ScaffoldDescriptorDTO:
    scaffoldId: int
    scaffoldName: str
    spacerLength: int
    scaffoldBordersBP: ScaffoldBordersBP

    @staticmethod
    def fromEntity(descriptor: Tuple[ScaffoldDescriptor, ScaffoldBordersBP]) -> 'ScaffoldDescriptorDTO':
        return ScaffoldDescriptorDTO(
            int(descriptor[0].scaffold_id),
            descriptor[0].scaffold_name,
            int(descriptor[0].spacer_length if descriptor[0].spacer_length is not None else 1000),
            ScaffoldBordersBPDTO.fromEntity(descriptor[1])
        )


@dataclass
class AssemblyInfo:
    contigDescriptors: List[Tuple[ContigDescriptor, ContigDirection]]
    scaffoldDescriptors: List[ScaffoldDescriptor]


@dataclass
class AssemblyInfoDTO:
    contigDescriptors: List[ContigDescriptorDTO]
    scaffoldDescriptors: List[ScaffoldDescriptorDTO]

    @staticmethod
    def fromEntity(assembly: AssemblyInfo) -> 'AssemblyInfoDTO':
        return AssemblyInfoDTO(
            [ContigDescriptorDTO.fromEntity(descriptor, dir)
             for descriptor, dir in assembly.contigDescriptors],
            [ScaffoldDescriptorDTO.fromEntity(descriptor)
             for descriptor in assembly.scaffoldDescriptors
             ]
        )


@dataclass
class GroupContigsIntoScaffoldRequest:
    start_bp: np.int64
    end_bp: np.int64
    name: Optional[str]
    spacer_length: Optional[int]


@dataclass
class GroupContigsIntoScaffoldRequestDTO:
    start_bp: int
    end_bp: int
    name: Optional[str]
    spacer_length: Optional[int]

    def __init__(self, request_json) -> None:
        self.start_bp: int = int(request_json['startBP'])
        self.end_bp: int = int(request_json['endBP'])
        self.name: Optional[str] = (
            request_json['scaffoldName'] if 'scaffoldName' in request_json.keys() else None)
        self.spacer_length: Optional[int] = int(
            request_json['spacerLength']) if 'spacerLength' in request_json.keys() else None

    def toEntity(self) -> GroupContigsIntoScaffoldRequest:
        return GroupContigsIntoScaffoldRequest(
            np.int64(self.start_bp),
            np.int64(self.end_bp),
            self.name if self.name != "" else None,
            self.spacer_length
        )


@dataclass
class UngroupContigsFromScaffoldRequest:
    start_bp: np.int64
    end_bp: np.int64


@dataclass
class UngroupContigsFromScaffoldRequestDTO:
    start_bp: int
    end_bp: int

    def __init__(self, request_json) -> None:
        self.start_bp: int = int(request_json['startBP'])
        self.end_bp: int = int(request_json['endBP'])

    def toEntity(self) -> UngroupContigsFromScaffoldRequest:
        return UngroupContigsFromScaffoldRequest(
            np.int64(self.start_bp),
            np.int64(self.end_bp),
        )


@dataclass
class ReverseSelectionRangeRequest:
    start_bp: np.int64
    end_bp: np.int64


@dataclass
class ReverseSelectionRangeRequestDTO:
    start_bp: int
    end_bp: int

    def __init__(self, request_json) -> None:
        self.start_bp: int = int(request_json['startBP'])
        self.end_bp: int = int(request_json['endBP'])

    def toEntity(self) -> ReverseSelectionRangeRequest:
        return ReverseSelectionRangeRequest(
            np.int64(self.start_bp),
            np.int64(self.end_bp),
        )


@dataclass
class MoveSelectionRangeRequest:
    start_bp: np.int64
    end_bp: np.int64
    target_start_bp: np.int64


@dataclass
class MoveSelectionRangeRequestDTO:
    start_bp: int
    end_bp: int
    target_start_bp: int

    def __init__(self, request_json) -> None:
        self.start_bp: int = int(request_json['startBP'])
        self.end_bp: int = int(request_json['endBP'])
        self.target_start_bp: int = int(request_json['targetStartBP'])

    def toEntity(self) -> MoveSelectionRangeRequest:
        return MoveSelectionRangeRequest(
            np.int64(self.start_bp),
            np.int64(self.end_bp),
            np.int64(self.target_start_bp)
        )
        
        
@dataclass
class SplitContigRequest:
    split_px: np.int64
    bp_resolution: np.int64


@dataclass
class SplitContigRequestDTO:
    split_px: int
    bp_resolution: int

    def __init__(self, request_json) -> None:
        self.split_px: int = int(request_json['splitPx'])
        self.bp_resolution: int = int(request_json['bpResolution'])

    def toEntity(self) -> SplitContigRequest:
        return SplitContigRequest(
            np.int64(self.split_px),
            np.int64(self.bp_resolution),
        )


@dataclass
class GetFastaForSelectionRequest:
    from_bp_x: np.int64
    from_bp_y: np.int64
    to_bp_x: np.int64
    to_bp_y: np.int64


@dataclass
class GetFastaForSelectionRequestDTO:
    from_bp_x: int
    from_bp_y: int
    to_bp_x:   int
    to_bp_y:   int

    def __init__(self, request_json) -> None:
        self.from_bp_x: int = int(request_json['fromBpX'])
        self.from_bp_y: int = int(request_json['fromBpY'])
        self.to_bp_x: int = int(request_json['toBpX'])
        self.to_bp_y: int = int(request_json['toBpY'])

    def toEntity(self) -> GetFastaForSelectionRequest:
        return GetFastaForSelectionRequest(
            np.int64(self.from_bp_x),
            np.int64(self.from_bp_y),
            np.int64(self.to_bp_x),
            np.int64(self.to_bp_y),
        )


@dataclass
class OpenFileResponse:
    status: str
    dtype: np.dtype
    resolutions: List[np.int64]
    pixelResolutions: List[np.float64]
    tileSize: int
    assemblyInfo: AssemblyInfo
    matrixSizesBins: List[np.int64]


@dataclass
class OpenFileResponseDTO:
    status: str
    dtype: str
    resolutions: List[int]
    pixelResolutions: List[float]
    tileSize: int
    assemblyInfo: AssemblyInfoDTO
    matrixSizesBins: List[int]

    @staticmethod
    def fromEntity(response: OpenFileResponse) -> 'OpenFileResponseDTO':
        return OpenFileResponseDTO(
            str(response.status),
            str(response.dtype),
            [int(r) for r in response.resolutions],
            [float(pr) for pr in response.pixelResolutions],
            int(response.tileSize),
            AssemblyInfoDTO.fromEntity(response.assemblyInfo),
            [int(sz) for sz in response.matrixSizesBins]
        )


@dataclass
class NormalizationSettings:
    preLogBase: np.float64
    preLogLnBase: np.float64
    postLogBase: np.float64
    postLogLnBase: np.float64
    applyCoolerWeights: bool

    def preLogEnabled(self) -> bool:
        return self.preLogBase > np.int64(1.0)

    def postLogEnabled(self) -> bool:
        return self.postLogBase > np.int64(1.0)

    def coolerBalanceEnabled(self) -> bool:
        return self.applyCoolerWeights


@dataclass
class NormalizationSettingsDTO:
    def __init__(self, request_json) -> None:
        self.preLogBase = request_json['preLogBase']
        self.postLogBase = request_json['postLogBase']
        self.applyCoolerWeights = request_json['applyCoolerWeights']

    def toEntity(self) -> NormalizationSettings:
        return NormalizationSettings(
            np.float64(self.preLogBase),
            np.log(np.float64(self.preLogBase)) if np.float64(
                self.preLogBase) > 1.0 else 1.0,
            np.float64(self.postLogBase),
            np.log(np.float64(self.postLogBase)) if np.float64(
                self.postLogBase) > 1.0 else 1.0,
            bool(self.applyCoolerWeights)
        )


@dataclass
class ContrastRangeSettings:
    lowerSignalBound: np.float64
    upperSignalBound: np.float64


@dataclass
class ContrastRangeSettingsDTO:
    def __init__(self, request_json) -> None:
        self.lowerSignalBound = request_json['lowerSignalBound']
        self.upperSignalBound = request_json['upperSignalBound']

    def toEntity(self) -> ContrastRangeSettings:
        return ContrastRangeSettings(
            np.float64(self.lowerSignalBound),
            np.float64(self.upperSignalBound)
        )
