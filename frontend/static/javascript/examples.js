const codeSnippets = {
    python: `import requests\nresponse = requests.get('https://api.example.com')\nprint(response.json())`,
    javascript: `fetch('https://api.example.com')\n  .then(response => response.json())\n  .then(data => console.log(data));`,
    java: `import java.net.*;\nimport java.io.*;\npublic class Main {\n  public static void main(String[] args) throws Exception {\n    URL url = new URL("https://api.example.com");\n    BufferedReader in = new BufferedReader(new InputStreamReader(url.openStream()));\n    String inputLine;\n    while ((inputLine = in.readLine()) != null)\n      System.out.println(inputLine);\n    in.close();\n  }\n}`,
    ruby: `require 'net/http'\nuri = URI('https://api.example.com')\nresponse = Net::HTTP.get(uri)\nputs response`,
    r: `library(httr)\nresponse <- GET('https://api.example.com')\ncontent(response, "text")`,
    perl: `use LWP::Simple;\nmy $content = get("https://api.example.com");\nprint $content;`
};

function showCode(language) {
    const codeElement = document.getElementById("code-snippet");
    codeElement.textContent = codeSnippets[language];

    // Highlight the active button
    const buttons = document.querySelectorAll(".btn-api");
    buttons.forEach((btn) => {
        btn.classList.remove("active");
        if (btn.innerText.toLowerCase() === language) {
            btn.classList.add("active");
        }
    });
}
