# Color palettes

get_pals <- function(pal) {
    if (!(pal %in% 1:9)) {
        stop("No palette matching supplied argument.")
    }
    switch(
        pal,
        c(
            "#44B3C2",
            "#F1A94E",
            "#E45641",
            "#5D4C46",
            "#7B8D8E",
            "#F2EDD8"
        ),
        c("#f4f6af", "#e3c2bd", "#d28eca", "#9b7cc3", "#656abb"),
        c("#011f4b", "#03396b", "#005b96", "#6497b1", "#b3cde0"),
        c(
            "#987543",
            "#E4C89F",
            "#BB9866",
            "#7D5A26",
            "#57390D",
            "#988443",
            "#E4D49F",
            "#BBA766",
            "#7D6926",
            "#57460D",
            "#3A376A",
            "#78769F",
            "#525082",
            "#262357",
            "#13103D",
            "#304A63",
            "#6A8094",
            "#466079",
            "#1D3751",
            "#0C2339"
        ),
        c(
            "#E58300",
            "#ADD8C6",
            "#DDAABB",
            "#664400",
            "#4C453D",
            "#002211",
            "#101111",
            "#881144",
            "#114488",
            "#775566",
            "#3F7F64",
            "#2B2B2B"
        ),
        c(
            "#ffffcc",
            "#c7e9b4",
            "#7fcdbb",
            "#41b6c4",
            "#2c7fb8",
            "#253494"
        ),
        c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6"),
        c(
            "#762a83",
            "#af8dc3",
            "#e7d4e8",
            "#d9f0d3",
            "#7fbf7b",
            "#1b7837"
        ),
        c(
            "#8c510a",
            "#d8b365",
            "#f6e8c3",
            "#c7eae5",
            "#5ab4ac",
            "#01665e"
        )
    )
}

display_pal <- function(pal) {
    barplot(rep(1, length(pal)), col = pal, axes = FALSE)
}
